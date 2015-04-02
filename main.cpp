#include <cmath>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

int usage() {
    cout << "\n\
Usage:   cross-sampler [options] <in.bam>\n\
Options: -t FILE    table of current and target read depth for each windows, required\n\
         -o FILE    output BAM\n\
         -u         output uncompressed BAMs\n\
         -s INT     seed for random shuffling\n\
         -f STR     required read flag\n\
         -F STR     filtered read flag\n\
         -?         print this message\n\
\n";
    return 1;
}

int read_region_depth(const char* fn_tgt, BamReader& reader, vector<BamRegion>& regions, vector<unsigned int>& src_depths, vector<unsigned int>& tgt_depths) {
    ifstream ifs;
    string line, chr;
    int chr_id, begin, end, src_dp, tgt_dp;

    ifs.open(fn_tgt);
    if (!ifs.is_open()) {
        cerr << "ERROR: failed to open [" << fn_tgt << "]\n";
        return 1;
    }
    while (getline(ifs, line)) {
        stringstream linestream(line);
        linestream >> chr >> begin >> end >> src_dp >> tgt_dp;
        chr_id = reader.GetReferenceID(chr);
        if (chr_id == -1) {
            cerr << "WARNING: sequence [" << chr << "] not found in header, region [" << chr << ':' << begin << '-' << end << "] skipped\n";
        }
        else {
            regions.push_back(BamRegion(chr_id, begin, chr_id, end));
            src_depths.push_back(src_dp);
            tgt_depths.push_back(tgt_dp);
        }
    }
    cerr << "INFO: [" << regions.size() << "] regions read\n";

    return 0;
}

int cigar_string(vector<CigarOp>& cigar_data, char* buf, int buf_size) {
    char* p=buf;
    int n=0;
    for (vector<CigarOp>::iterator it=cigar_data.begin(); it!=cigar_data.end(); ++it) {
        n = floor(log10(it->Length))+2;
        if (p+n < buf+buf_size-1) sprintf(p, "%d%c", it->Length, it->Type);
        else return p+n-buf;
        p += n;
    }
    *p = '\0';
    return 0;
}

void print_alignment(BamAlignment& aln, const RefVector& refseq, FILE* file=stdout) {
    char cigar[100];
    if (cigar_string(aln.CigarData, cigar, 100) != 0) {
        cerr << "WARNING: failed to compile cigar string for [" << aln.Name.c_str() << "]\n";
        return;
    }
    string mate_chr = refseq[aln.MateRefID].RefName;
    if (mate_chr.compare(refseq[aln.RefID].RefName) == 0) mate_chr = "=";
    fprintf(file, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
            aln.Name.c_str(),
            aln.AlignmentFlag,
            refseq[aln.RefID].RefName.c_str(),
            aln.Position,
            aln.MapQuality,
            cigar,
            mate_chr.c_str(),
            aln.MatePosition,
            aln.InsertSize,
            aln.QueryBases.c_str(),
            aln.Qualities.c_str()
            );
}

int get_overlap(BamAlignment& aln, BamRegion& region) {
    int ovlp_beg = max(aln.Position, region.LeftPosition);
    int ovlp_end = min(aln.GetEndPosition(), region.RightPosition);
    return max(0, ovlp_end-ovlp_beg);
}

int main(const int argc, char* const argv[]) {
    int c, min_mapQ=0, seed=chrono::system_clock::now().time_since_epoch().count();
    unsigned int flag_on=0, flag_off=0;
    string fn_tgt, fn_in, fn_out="", out_format="b";

    while ((c = getopt(argc, argv, "SbBcCt:h1Ho:q:f:F:ul:r:?T:R:L:s:@:m:x:U:")) >= 0) {
        switch (c) {
            case 's': seed = atoi(optarg); break;
            case 'm': break;
            case 'c': break;
            case 'S': break;
            case 'b': break;
            case 'C': break;
            case 'h': break;
            case 'H': break;
            case 'o': fn_out = optarg; break;
            case 'U': break;
            case 'f': flag_on |= strtol(optarg, 0, 0); break;
            case 'F': flag_off |= strtol(optarg, 0, 0); break;
            case 'q': min_mapQ = atoi(optarg); break;
            case 'u': out_format = "u"; break;
            case '1': break;
            case 'l': break;
            case 'r': break;
            case 't': fn_tgt = optarg; break;
            case 'R': break;
            case '?': return usage();
            case 'T': break;
            case 'B': break;
            case '@': break;
            case 'x': break;
            default: return usage();
        }
    }
    if (fn_tgt.compare("") == 0) return usage();
    if (argc == optind) return usage();
    fn_in = argv[optind];

    BamReader reader;
    if (!reader.Open(fn_in)) {
        cerr << "ERROR: cannot open [" << fn_in << "] for reading\n";
        return 1;
    }
    if (!reader.LocateIndex()) {
        cerr << "ERROR: cannot find BAM index for [" << fn_in << "]\n";
        return 1;
    }

    const SamHeader header = reader.GetHeader();
    if (header.SortOrder.compare("coordinate") != 0) {
        cerr << "ERROR: [" << fn_in << "] not sorted by coordinate\n";
        return 1;
    }
    const RefVector refseq = reader.GetReferenceData();

    vector<BamRegion> regions;
    vector<unsigned int> src_depths, tgt_depths;
    if (read_region_depth(fn_tgt.c_str(), reader, regions, src_depths, tgt_depths) != 0) return 1;

    BamWriter writer;
    if (!writer.Open(fn_out, header, refseq)) {
        cerr << "ERROR: cannot open [" << fn_out << "] for writing\n";
        return 1;
    }

    BamAlignment aln;
    vector<BamAlignment> reads;
    vector<string> paired, unpaired;
    unordered_map<int, int> kept;
    unordered_map<string, unsigned int> seen, sampled;
    unordered_map<string, vector<int> > pool;
    for (size_t i=0; i<regions.size(); ++i) {
        reads.clear();
        paired.clear();
        unpaired.clear();
        kept.clear();
        pool.clear();

        char region_string[256];
        sprintf(region_string, "%s:%d-%d", refseq[regions[i].LeftRefID].RefName.c_str(), regions[i].LeftPosition, regions[i].RightPosition);
        if (!reader.SetRegion(regions[i])) {
            cerr << "WARNING: failed to locate [" << region_string << "]\n";
            //cerr << "WARNING: failed to locate [" << refseq[regions[i].LeftRefID].RefName << ':' << regions[i].LeftPosition << '-' << regions[i].RightPosition << "]\n";
            continue;
        }
        while (reader.GetNextAlignment(aln)) {
            if ((aln.AlignmentFlag & flag_on) == flag_on && !(aln.AlignmentFlag & flag_off) && aln.MapQuality >= min_mapQ)
                reads.push_back(aln);
        }
        if (reads.size() == 0) continue;

        unsigned int depth = 0;
        for (size_t k=0; k<reads.size(); ++k) {
            aln = reads[k];
            string rn = aln.Name;
            if (seen.find(rn) != seen.end()) { // if seen in previous regions
                if (sampled.find(rn) != sampled.end()) { // if self or mate sampled before, sample it
                    if (sampled[rn] != aln.AlignmentFlag) kept[k] = 1; // if mate sampled before, keep it
                    depth += get_overlap(aln, regions[i]);
                }
                if (seen[rn] != aln.AlignmentFlag) seen[rn] = aln.AlignmentFlag;
            }
            else { // if not seen in previous regions
                pool[rn].push_back(k);
            }
            if (depth > tgt_depths[i]) break;
        }
        if (depth < tgt_depths[i]) {
            for (auto it=pool.begin(); it!=pool.end(); ++it) {
                if (it->second.size()>1)
                    paired.push_back(it->first);
                else
                    unpaired.push_back(it->first);
            }
            shuffle(paired.begin(), paired.end(), default_random_engine(seed));
            shuffle(unpaired.begin(), unpaired.end(), default_random_engine(seed));
            int n1=paired.size(), n2=unpaired.size(), k1, k2, k3;
            while (depth < tgt_depths[i] && n1+n2 > 0) {
                if (n1>0) {
                    k1 = pool[paired[--n1]][0];
                    k2 = pool[paired[n1]][1];
                    depth += get_overlap(reads[k1], regions[i]);
                    depth += get_overlap(reads[k2], regions[i]);
                    kept[k1] = 1; kept[k2] = 1;
                    continue;
                }
                if (n2>0) {
                    k3 = pool[unpaired[--n2]][0];
                    depth += get_overlap(reads[k3], regions[i]);
                    kept[k3] = 1;
                    continue;
                }
            }
        }
        for (auto it=pool.begin(); it!=pool.end(); ++it) {
            string rn = it->first;
            seen[rn] = reads[pool[rn].back()].AlignmentFlag;
        }
        for (auto it=kept.begin(); it!=kept.end(); ++it) {
            int k = it->first;
            string rn = reads[k].Name;
            sampled[rn] = reads[k].AlignmentFlag;
            writer.SaveAlignment(reads[k]);
        }
        cerr << "INFO: target=[" << tgt_depths[i] << "], actual=[" << depth << "], N(reads)=[" << reads.size() << "], N(kept)=[" << kept.size() << "] at [" << region_string << "]\n";
    }

    reader.Close();

    return 0;
}
