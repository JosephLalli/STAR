#include "PersonalizedInput.h"

#include "embedded_consensus/ConsensusEmbedded.h"
#include "htslib/htslib/synced_bcf_reader.h"
#include "htslib/htslib/vcf.h"

#include <cctype>
#include <climits>
#include <cstdlib>
#include <exception>
#include <memory>
#include <stdexcept>

namespace {
string uppercase(const char *bases)
{
    if (bases == nullptr) throw runtime_error("malformed VCF/BCF allele");
    string result(bases);
    for (char &base : result)
        base = static_cast<char>(toupper(static_cast<unsigned char>(base)));
    return result;
}

bool concrete(const string &bases)
{
    if (bases.empty()) return false;
    for (char base : bases)
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N')
            return false;
    return true;
}

bool symbolic(const string &bases)
{
    return bases.find_first_of("<>[]*") != string::npos;
}

string readerError(const bcf_srs_t *reader)
{
    const char *message = reader == nullptr ? nullptr : bcf_sr_strerror(reader->errnum);
    return message == nullptr ? string("unknown HTSlib error") : string(message);
}

string regionName(const string &contig)
{
    return contig.find_first_of(",:-{}") == string::npos ? contig : "{" + contig + "}";
}

} // namespace

// Own one successfully opened reader. The failed-add workaround below remains tied to the pinned HTS ABI.
struct PersonalizedInput::VcfReader {
    bcf_srs_t *reader;
    ~VcfReader() { bcf_sr_destroy(reader); }
};

void PersonalizedInput::closeVcfReader()
{
    if (vcfReader != nullptr) {
        delete vcfReader;
        vcfReader = nullptr;
    }
}

PersonalizedInput::VcfReader *PersonalizedInput::openVcfReader() const
{
    bcf_srs_t *reader = bcf_sr_init();
    if (reader == nullptr) throw runtime_error("cannot initialize HTSlib VCF/BCF reader");
    if (bcf_sr_set_opt(reader, BCF_SR_REQUIRE_IDX) != 0) {
        bcf_sr_destroy(reader);
        throw runtime_error("cannot require indexed personalization VCF/BCF input");
    }
    if (!bcf_sr_add_reader(reader, vcfPath.c_str())) {
        const string error = readerError(reader);
        // HTSlib's failed add owns an incomplete reader slot. Avoid its
        // destructor dereferencing that incomplete slot on this fatal path.
        reader->nreaders = 0;
        bcf_sr_destroy(reader);
        throw runtime_error("personalization VCF/BCF requires indexed BGZF VCF.gz or BCF: " + error);
    }
    if (reader->errnum != 0) {
        const string error = readerError(reader);
        bcf_sr_destroy(reader);
        throw runtime_error("personalization VCF/BCF input may be truncated: " + error);
    }

    bcf_hdr_t *header = bcf_sr_get_header(reader, 0);
    if (header == nullptr) {
        bcf_sr_destroy(reader);
        throw runtime_error("cannot read personalization VCF/BCF header");
    }
    const int selected = bcf_hdr_set_samples(header, vcfSample.c_str(), 0);
    if (selected != 0 || bcf_hdr_nsamples(header) != 1) {
        bcf_sr_destroy(reader);
        throw runtime_error("personalization sample absent from VCF/BCF: " + vcfSample);
    }
    try { return new VcfReader{reader}; }
    catch (...) { bcf_sr_destroy(reader); throw; }
}

vector<PersonalizedInput::SelectedCall> PersonalizedInput::selectCalls(
    VcfReader *owner, const string &contig, const string &reference,
    map<string,uint64> &counts) const
{
    if (owner == nullptr) throw runtime_error("personalization VCF/BCF reader is not available");
    bcf_srs_t *reader = owner->reader;
    if (reference.empty() || reference.size() > static_cast<size_t>(INT_MAX))
        throw runtime_error("reference contig must have a positive signed-32-bit length: " + contig);

    bcf_hdr_t *header = bcf_sr_get_header(reader, 0);
    if (header == nullptr || bcf_hdr_nsamples(header) != 1)
        throw runtime_error("personalization reader lost its selected donor sample");

    vector<SelectedCall> selectedCalls;
    if (bcf_hdr_name2id(header, contig.c_str()) < 0) return selectedCalls;

    const string region = regionName(contig);
    if (bcf_sr_set_regions(reader, region.c_str(), 0) != 0)
        throw runtime_error("cannot seek indexed personalization VCF/BCF contig: " + contig);
    if (reader->errnum != 0)
        throw runtime_error("cannot seek indexed personalization VCF/BCF contig " + contig + ": " + readerError(reader));

    // bcftools consensus freezes each applied REF span through its last base and skips a later call starting
    // inside it, except an indel anchored on that last base after a non-insertion (D091).
    array<int,2> frozen = {{-1, -1}};
    array<bool,2> previousInsertion = {{false, false}};
    int previousPosition = -1;
    while (bcf_sr_next_line(reader)) {
        if (!bcf_sr_has_line(reader, 0))
            throw runtime_error("personalization reader returned no record for its only input");
        bcf1_t *record = bcf_sr_get_line(reader, 0);
        if (record == nullptr) throw runtime_error("personalization reader returned a null record");
        bcf_unpack(record, BCF_UN_STR);

        ++counts["records_scanned"];
        if (record->pos < 0 || record->pos > INT_MAX)
            throw runtime_error("personalization VCF/BCF coordinate is outside signed-32-bit range on " + contig);
        const int position = static_cast<int>(record->pos);
        if (position < previousPosition)
            throw runtime_error("personalization VCF/BCF records are not coordinate sorted on " + contig);
        previousPosition = position;

        int32_t *rawGenotypes = nullptr;
        int allocated = 0;
        const int genotypeCount = bcf_get_genotypes(header, record, &rawGenotypes, &allocated);
        std::unique_ptr<int32_t, decltype(&free)> genotypeOwner(rawGenotypes, &free);
        vector<int32_t> genotype;
        if (genotypeCount > 0) {
            bool vectorEnded = false;
            for (int index = 0; index < genotypeCount; ++index) {
                const int32_t encoded = rawGenotypes[index];
                if (encoded == bcf_int32_vector_end) {
                    vectorEnded = true;
                    continue;
                }
                if (vectorEnded)
                    throw runtime_error("malformed GT vector has values after vector-end on " + contig);
                genotype.push_back(encoded);
            }
        }

        if (genotype.empty()) {
            ++counts["no_carried_alt_records"];
            continue;
        }
        bool carriesAlternative = false;
        for (int32_t encoded : genotype)
            if (!bcf_gt_is_missing(encoded) && bcf_gt_allele(encoded) > 0)
                carriesAlternative = true;
        if (!carriesAlternative) {
            ++counts["no_carried_alt_records"];
            continue;
        }

        ++counts["nonreference_records"];
        ++counts["ledger_rows"];
        if (genotype.size() != 1 && genotype.size() != 2)
            throw runtime_error("unsupported genotype ploidy on " + contig);
        // As bcftools consensus -H 1 / -H 2 select: each haplotype takes the GT allele in its slot, in
        // written order whether phased or not. A missing allele, the REF allele, or a '*'/symbolic
        // allele leaves only that haplotype unchanged: bcftools writes an orphan '*' literally, and STAR
        // would read it as N (D026, D090).
        array<int,2> alleles;
        array<string,2> alternates;
        array<bool,2> concreteAlt = {{false, false}};
        for (size_t haplotype = 0; haplotype < 2; ++haplotype) {
            const int32_t encoded = genotype[genotype.size() == 1 ? 0 : haplotype];
            alleles[haplotype] = bcf_gt_is_missing(encoded) ? -1 : bcf_gt_allele(encoded);
            if (alleles[haplotype] >= record->n_allele)
                throw runtime_error("GT allele index beyond ALT list on " + contig);
            if (alleles[haplotype] > 0) {
                alternates[haplotype] = uppercase(record->d.allele[alleles[haplotype]]);
                concreteAlt[haplotype] = !symbolic(alternates[haplotype]);
            }
        }
        if (genotype.size() == 2 && alleles[0] >= 0 && alleles[1] >= 0 && alleles[0] != alleles[1]
                && !bcf_gt_is_phased(genotype[1]))
            ++counts["unphased_het_records_in_gt_order"];
        const string ref = uppercase(record->d.allele[0]);
        if (ref.empty()) throw runtime_error("empty REF allele on " + contig);
        if ((concreteAlt[0] || concreteAlt[1])
                && static_cast<uint64>(position) + ref.size() > reference.size()) {
            ++counts["h1_ref_span_past_contig_end_skipped"];
            ++counts["h2_ref_span_past_contig_end_skipped"];
            continue;
        }
        if ((concreteAlt[0] || concreteAlt[1])
                && reference.compare(static_cast<size_t>(position), ref.size(), ref) != 0)
            throw runtime_error("REF mismatch in personalization VCF/BCF at " + contig + ":" + to_string(position + 1));

        SelectedCall selected;
        selected.position = position;
        selected.reference = ref;
        selected.applied = {{false, false}};
        selected.sharesAnchor = {{false, false}};
        for (size_t haplotype = 0; haplotype < alleles.size(); ++haplotype) {
            const int allele = alleles[haplotype];
            const string prefix = string("h") + to_string(haplotype + 1);
            if (allele < 0) {
                ++counts[prefix + "_missing_allele_skipped"];
                continue;
            }
            if (allele == 0) {
                ++counts[prefix + "_reference_allele"];
                continue;
            }
            if (!concreteAlt[haplotype]) {
                ++counts[prefix + "_star_or_symbolic_skipped"];
                continue;
            }
            const string &alternate = alternates[haplotype];
            bool sharesAnchor = false;
            if (position <= frozen[haplotype]) {
                sharesAnchor = position == frozen[haplotype] && !previousInsertion[haplotype]
                    && (bcf_get_variant_type(record, allele) & VCF_INDEL)
                    && record->d.allele[0][0] == record->d.allele[allele][0]
                    && record->d.var[allele].n != 0;
                if (!sharesAnchor) {
                    ++counts[prefix + "_overlap_skipped"];
                    continue;
                }
            }
            if (alternate == ref) {
                ++counts[string("h") + to_string(haplotype + 1) + "_noop_after_trim"];
                continue;
            }
            if (!concrete(ref) || !concrete(alternate))
                throw runtime_error("unsupported concrete sequence alphabet on " + contig + ":" + to_string(position + 1));
            selected.alternate[haplotype] = alternate;
            selected.applied[haplotype] = true;
            selected.sharesAnchor[haplotype] = sharesAnchor;
            frozen[haplotype] = position + static_cast<int>(ref.size()) - 1;
            previousInsertion[haplotype] = alternate.size() > ref.size();
            ++counts[prefix + "_applied"];
            if (sharesAnchor) ++counts[prefix + "_applied_shared_anchor"];
        }
        if (selected.applied[0] || selected.applied[1]) selectedCalls.push_back(std::move(selected));
    }
    if (reader->errnum != 0)
        throw runtime_error("failed while reading personalization VCF/BCF contig " + contig + ": " + readerError(reader));
    return selectedCalls;
}

PersonalizedInput::FastaResult PersonalizedInput::prepareContig(
    const string &source, const string &reference) const
{
    FastaResult result;
    try {
        unique_ptr<VcfReader> reader(openVcfReader());
        const vector<SelectedCall> calls = selectCalls(reader.get(), source, reference, result.counts);

        const int outputHaplotypes = haploid.count(source) == 0 ? 2 : 1;
        result.sequences.reserve(static_cast<size_t>(outputHaplotypes));
        if (!skipAnnotation.count(source))
            result.haplotypes.reserve(static_cast<size_t>(outputHaplotypes));
        for (int haplotype = 0; haplotype < outputHaplotypes; ++haplotype) {
            const string suffix = outputHaplotypes == 1 ? string() : (haplotype == 0 ? "_L" : "_R");
            const string outputName = source + suffix;
            string sequence;
            sequence.reserve(reference.size());
            ConsensusEmbedded *engine = ce_create(source.c_str(), static_cast<int>(reference.size()));
            if (engine == nullptr) throw runtime_error("cannot initialize embedded consensus for " + source);
            try {
                int cursor = 0;
                for (const SelectedCall &call : calls) {
                    if (!call.applied[haplotype]) continue;
                    // A shared-anchor call starts on the previous call's last REF base, already emitted by it.
                    const bool sharesAnchor = call.sharesAnchor[haplotype];
                    if (sharesAnchor ? call.position != cursor - 1 : call.position < cursor)
                        throw runtime_error("selected overlapping calls reached consensus engine on " + source);
                    if (!sharesAnchor)
                        sequence.append(reference.data() + cursor,
                                        static_cast<size_t>(call.position - cursor));
                    const char *edited = nullptr;
                    int editedLength = 0;
                    if (ce_apply(engine, call.position, call.reference.c_str(),
                                 call.alternate[haplotype].c_str(), &edited, &editedLength) != 0
                            || edited == nullptr || editedLength < 0)
                        throw runtime_error("embedded consensus failed on selected VCF/BCF call at "
                                            + source + ":" + to_string(call.position + 1));
                    sequence.append(edited, static_cast<size_t>(editedLength));
                    cursor = call.position + static_cast<int>(call.reference.size());
                }
                sequence.append(reference.data() + cursor,
                                reference.size() - static_cast<size_t>(cursor));
                ce_finish(engine);

                if (sequence.empty() || sequence.size() > static_cast<size_t>(INT_MAX))
                    throw runtime_error("personalized contig exceeds supported signed-32-bit length: " + source);
                if (!skipAnnotation.count(source)) {
                    Haplotype output;
                    output.name = outputName;
                    output.suffix = suffix;
                    output.referenceSize = static_cast<int>(reference.size());
                    output.querySize = static_cast<int>(sequence.size());
                    const int blockCount = ce_block_count(engine);
                    if (blockCount < 0)
                        throw runtime_error("embedded consensus did not finish its chain for " + source);
                    output.chain.reserve(static_cast<size_t>(blockCount));
                    int previousEnd = 0;
                    int previousQueryEnd = 0;
                    for (int block = 0; block < blockCount; ++block) {
                        ChainBlock mapped;
                        ce_block(engine, block, &mapped.referenceStart, &mapped.referenceEnd,
                                 &mapped.queryStart);
                        if (mapped.referenceStart < previousEnd
                                || mapped.referenceEnd <= mapped.referenceStart
                                || mapped.referenceEnd > output.referenceSize
                                || mapped.queryStart < previousQueryEnd
                                || static_cast<int64>(mapped.queryStart) + mapped.referenceEnd
                                       - mapped.referenceStart > output.querySize)
                            throw runtime_error("embedded consensus emitted an invalid forward chain for " + source);
                        output.chain.push_back(mapped);
                        previousEnd = mapped.referenceEnd;
                        previousQueryEnd = mapped.queryStart + mapped.referenceEnd - mapped.referenceStart;
                    }
                    result.haplotypes.push_back(std::move(output));
                }
            } catch (...) {
                ce_destroy(engine);
                throw;
            }
            ce_destroy(engine);
            result.sequences.push_back(std::move(sequence));
        }
    } catch (const std::exception &error) {
        result.error = error.what();
    } catch (...) {
        result.error = "unknown error while preparing personalized contig " + source;
    }
    return result;
}
