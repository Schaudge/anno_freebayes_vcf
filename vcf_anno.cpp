/* 
 *  Annotate the snp/indel variant vcf records from freebayes !
 *  Written By Schaudge King
 *  Date: 2017-01-10
 */
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include "log/gcbiLogCpp.h"
#include "anno/CallRefDatabase.h"
#include "anno/DenovoInvoke.h"
#include "anno/RefSnpTable.h"
#include "anno/FreqAnnoTable.h"

using namespace std;

// void writtenAllRecord(ofstream&, const string&, const vector< tuple<string,string> >&);
void mappingVarInfor(ifstream&, ofstream&, RefSnpTable &,  FreqAnnoTable &, log4cpp::Category &, const string command="/gcbi/analyze/integration/denovo.sh");
void mappingVarInfor(ifstream&, ofstream&, RefSnpTable &,  FreqAnnoTable &, log4cpp::Category &, sql::Connection * & , const string command="/gcbi/analyze/integration/denovo.sh");
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();		        
	}
	return str;
}
void insertToDatabase(sql::Connection * & mysqlConn, const std::vector< std::tuple<std::string, std::string> >& knownVar, const std::string& remainedFields) {
	std::string underChangedFields = ReplaceAll(remainedFields, std::string("\t"), std::string(","));
	std::string syntax = "REPLACE INTO denovo_brca (search_id,snp_id,chr,pos,ref,alt" + underChangedFields + ") VALUES (";
	sql::Statement * stmt = mysqlConn->createStatement();
	auto iter = knownVar.cbegin(), iter_end = knownVar.cend();
	for ( ; iter != iter_end; ++iter ) {
		std::string prefix = ReplaceAll(std::get<0>(*iter), std::string("\t"), std::string("\",\""));
		std::string postfix = ReplaceAll(std::get<1>(*iter), std::string("\t"), std::string("\",\""));
		//std::cout << "REPLACE Mysql ---> " << syntax << "\"" << prefix << postfix << "\");" << std::endl;
		stmt->execute(syntax + "\"" + prefix + postfix + "\");");
	}
}	

int main ( int argc , char * argv[] )
{
	int index; bool db_ref(0);
	string ref_anno{""}, vcf_file{""}, out_file{"1.vcf.anno"};
	string log_conf("/home/yuanshenran/cpp-workspace/vcfannotate/src/log/log4cpp.conf");
	while ( 0 <= (index = getopt(argc, argv, "r:v:l:o:d")) ) {
		switch (index) {
		case 'd':
			db_ref = 1;
			break;
		case 'r':
			ref_anno = optarg;
			break;
		case 'v':
			vcf_file = optarg;
			break;
		case 'o':
			out_file = optarg;
			break;
		case 'l':
			log_conf = optarg;
			break;
		default:
			fprintf(stderr, "Unrecognized option '-%c'.\n", index);
			return 1;
		}
	}
	if ( 1 == argc || 4 > argc || optind > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: vcf_annotation -r reference_anno -v vcf_file -l log_file_config\n");
		fprintf(stderr, "Program: vcf_annotation -d -v vcf_file -o vcf.anno");
		fprintf(stderr, "\n");
		return 1;
	}
	gcbiLogCpp logCate(log_conf);
	log4cpp::Category& Logger= logCate.getRefLog();
	Logger.info("Begin the vcf annotation process ... \t SAMPLE: " + vcf_file);

	// CallRefDatabase RefDatabase = CallRefDatabase("tcp://192.168.2.10", "devuser", "111111", "chgc-gcbi");
	CallRefDatabase RefDatabase = CallRefDatabase("tcp://127.0.0.1", "gcbi", "gCiB35G", "chgc-gcbi");
	ifstream vcf_stream(vcf_file);
	ofstream anno_file(out_file);
	// Loading the reference index database to memory(heap)!
	if ( db_ref ) {
		RefSnpTable recordDb = RefSnpTable(RefDatabase.getDbConnection());
		FreqAnnoTable freqBrcaDb = FreqAnnoTable(RefDatabase.getDbConnection());
		//freqBrcaDb.display_all();
		mappingVarInfor(vcf_stream, anno_file, recordDb, freqBrcaDb, Logger, RefDatabase.getDbConnection());
	} else {
		ifstream ref_snp_file(ref_anno);
		RefSnpTable recordDb = RefSnpTable(ref_snp_file);
		FreqAnnoTable freqBrcaDb = FreqAnnoTable(RefDatabase.getDbConnection());
		mappingVarInfor(vcf_stream, anno_file, recordDb, freqBrcaDb, Logger);
		ref_snp_file.close();
	}
	vcf_stream.close();
	anno_file.close();
	Logger.info("Now, the whole vcf annotation was completed.");
	return 0;
}


void writtenAllRecord (ofstream& outFile, const string& vcfInfor, FreqAnnoTable & freqDb, const vector< tuple<string,string> >& resultRecordSet ) 
{
	auto iter = resultRecordSet.cbegin(), iter_end = resultRecordSet.cend();
	string prefix, search_id;
	for ( ; iter != iter_end; ++iter ) {
		prefix = get<0>(*iter);
		search_id = prefix.substr(0, prefix.find_first_of("\t"));
		outFile << get<0>(*iter) << "\t" << vcfInfor << get<1>(*iter)
				   << ( (*freqDb.getFreqTable()).count(search_id) > 0 ? (*freqDb.getFreqTable()).at(search_id) : "\t\t\t\t\t" ) << "\n";
	}
}

void mappingVarInfor(ifstream & vcf_stream, ofstream & anno_file, RefSnpTable & recordDb, FreqAnnoTable & freqDb, log4cpp::Category & Logger, const string command)
{

	string line, field;
	string chr, pos, ref, alt, quality;
	vector<tuple<string,string> > knownVar;
	anno_file << "search_id\tsnp_id\tchr\tpos\tref\talt\tquality\tformat\tinfor\tallels_count" << recordDb.getOtherAnnoFields() << freqDb.getOtherFields() << "\n";
	while ( getline(vcf_stream,line) ) {
		if ( line[0] != '#' ) {
			auto i=0, j=0;
			istringstream record(line);
			record >> chr; record >> pos; record >> field;
			record >> ref; record >> alt; record >> quality;
			while (i<10) { i++; record >> field; }
			if ( chr[0] == 'c' ) chr.erase(0,3);
			if ( chr.size() < 2 ) chr = "0" + chr;

			string genotype, depth, RO, AO;
			i = field.find(":", j);
			genotype = field.substr(j, i);
			j = field.find(":", i+1);
			depth = field.substr(i+1,j-i-1);
			i = field.find(":",j+1); j = field.find(":", i+1);
			RO = field.substr(i+1, j-i-1);
			i = field.find(":",j+1); j = field.find(":", i+1);
			if ( stod(quality) <  1.0 || stoi(depth) < 10 ) continue;
			if ( genotype != "1/2" ) {
				AO = field.substr(i+1, j-i-1);
				string infor = quality+"\tGT:DP\t"+genotype+":"+depth+"\t"+RO+":"+AO;
				knownVar = getRelatedRecords(chr,pos,ref,alt,recordDb,true);
				if ( knownVar.size() >= 1 ) {
					writtenAllRecord(anno_file, infor, freqDb, knownVar);
				} else { // denovo variation site
					string novoSite = chr+"\t"+pos+"\t"+ref+"\t"+alt;
					DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
					if ( 0 == unknownSites.invoke_command(command) ) {
						knownVar = unknownSites.readResult();
						writtenAllRecord(anno_file, infor, freqDb, knownVar);
					} else {
						Logger.error("Unannotated denovo site: " + chr + " " + pos + " " + ref + " " + alt);
					}
				}
			} else {
				string ALL_AO = field.substr(i+1, j-i-1);
				i = ALL_AO.find(",");
				string AO1 = ALL_AO.substr(0,i), AO2 = ALL_AO.substr(i+1);
				string infor1 = quality+"\tGT:DP\t"+genotype+":"+depth+"\t"+RO+":"+AO1;
				string infor2 = quality+"\tGT:DP\t"+genotype+":"+depth+"\t"+RO+":"+AO2;
				// cout << infor1 << " -------- " << infor2 << endl;
				i = alt.find(",");
				knownVar = getRelatedRecords(chr,pos,ref,alt.substr(0,i),recordDb,false);
				if ( knownVar.size() >= 1 ) {
					writtenAllRecord(anno_file, infor1, freqDb, knownVar);
				} else { // denovo variation site
					string novoSite = chr+"\t"+pos+"\t"+ref+"\t"+alt.substr(0,i);
					DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
					if ( 0 == unknownSites.invoke_command(command) ) {
						knownVar = unknownSites.readResult();
						writtenAllRecord(anno_file, infor1, freqDb, knownVar);
					} else {
						Logger.error("Unannotated denovo site: " + chr + " " + pos + " " + ref + " " + alt);
					}
				}
				knownVar = getRelatedRecords(chr,pos,ref,alt.substr(i+1),recordDb,true);
				if ( knownVar.size() >= 1 ) {
					writtenAllRecord(anno_file, infor2, freqDb, knownVar);
				} else { // denovo variation site
					string novoSite = chr+"\t"+pos+"\t"+ref+"\t"+alt.substr(i+1);
					DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
					if ( 0 == unknownSites.invoke_command(command) ) {
						knownVar = unknownSites.readResult();
						writtenAllRecord(anno_file, infor2, freqDb, knownVar);
					} else {
						Logger.error("Unannotated denovo site: " + chr + " " + pos + " " + ref + " " + alt);
					}
				}
			}
		}
	}
}

void mappingVarInfor(ifstream & vcf_stream, ofstream & anno_file, RefSnpTable & recordDb, FreqAnnoTable & freqDb, log4cpp::Category & Logger, sql::Connection * & mysqlConn, const string command)
{

	string line, field;
	string chr, pos, ref, alt, quality;
	vector<tuple<string,string> > knownVar;
	anno_file << "search_id\tsnp_id\tchr\tpos\tref\talt\tquality\tformat\tinfor\tallels_count" << recordDb.getOtherAnnoFields() << freqDb.getOtherFields() << "\n";
	while ( getline(vcf_stream,line) ) {
		if ( line[0] != '#' ) {
			auto i=0, j=0;
			istringstream record(line);
			record >> chr; record >> pos; record >> field;
			record >> ref; record >> alt; record >> quality;
			while (i<10) { i++; record >> field; }
			if ( chr[0] == 'c' ) chr.erase(0,3);
			if ( chr.size() < 2 ) chr = "0" + chr;

			string genotype, depth, RO, AO;
			i = field.find(":", j);
			genotype = field.substr(j, i);
			j = field.find(":", i+1);
			depth = field.substr(i+1,j-i-1);
			i = field.find(":",j+1); j = field.find(":", i+1);
			RO = field.substr(i+1, j-i-1);
			i = field.find(":",j+1); j = field.find(":", i+1);
			if ( stod(quality) <  1.0 || stoi(depth) < 10 ) continue;
			if ( genotype != "1/2" ) {
				AO = field.substr(i+1, j-i-1);
				string infor = quality+"\tGT:DP\t"+genotype+":"+depth+"\t"+RO+":"+AO;
				std::cout << chr << " " << pos << " " << ref << "---" << infor << std::endl;
				knownVar = getRelatedRecords(chr,pos,ref,alt,recordDb,true);
				if ( knownVar.size() >= 1 ) {
					writtenAllRecord(anno_file, infor, freqDb, knownVar);
				} else { // denovo variation site
					string novoSite = chr+"\t"+pos+"\t"+ref+"\t"+alt;
					DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
					if ( 0 == unknownSites.invoke_command(command) ) {
						knownVar = unknownSites.readResult();
						writtenAllRecord(anno_file, infor, freqDb, knownVar);
						insertToDatabase(mysqlConn, knownVar, recordDb.getOtherAnnoFields());
					} else {
						Logger.error("Unannotated denovo site: " + chr + " " + pos + " " + ref + " " + alt);
					}
				}
			} else {
				string ALL_AO = field.substr(i+1, j-i-1);
				i = ALL_AO.find(",");
				string AO1 = ALL_AO.substr(0,i), AO2 = ALL_AO.substr(i+1);
				string infor1 = quality+"\tGT:DP\t0|1:"+depth+"\t"+RO+":"+AO1;
				string infor2 = quality+"\tGT:DP\t0|1:"+depth+"\t"+RO+":"+AO2;
				string infor = quality+"\tGT:DP\t1|1:"+depth+"\t"+RO+":"+std::to_string(stoul(AO1) + stoul(AO2));
				vector< tuple<string, string, string, string> > alterAllelesVector = { std::make_tuple("","","","") };
				map< string, vector< tuple<string, string, string, string> > > posVarMap;
				i = alt.find(",");
				string varPos{""};
				// two alternative alleles aligned to reference sequence !!!
				string ref_ = "", alt_ = "";
				string ref1_ = ref, ref2_ = ref, alt1_ = alt.substr(0,i), alt2_ = alt.substr(i+1);
				unsigned long max = ref1_.size() >= alt1_.size() ? ref1_.size() : alt1_.size();
				unsigned long min = ref1_.size() <= alt1_.size() ? ref1_.size() : alt1_.size();
				for ( auto j = max-1; j >= 0; j-- ) {
					if ( ref1_.size() == alt1_.size() ) {
						break;
					} else if ( ref1_.size() < alt1_.size() ) {
						auto coor = min - ( max - j );
						if ( ref1_[coor] != alt1_[j] ) {
							ref1_.insert(coor+1, "-");
							min++;
						}
					} else {
						auto coor = min - ( max - j );
						if ( alt1_[coor] != ref1_[j] ) {
							alt1_.insert(coor+1, "-");
							min++;
						}
					}
				}
				// Generate the separated variations !
				unsigned long ii = 0, jj = 1, inpos = 0;
				while ( ii < max ) {
					if ( ref1_[ii] != alt1_[ii] && ii+1 < max && ref1_.substr(ii+1,1) == string{"-"} ) {         // multiple sites insert variation
						jj = ii + 1;
						while ( jj+1 < max && ref1_.substr(jj+1,1) == string{"-"} ) jj++;
						if ( ref1_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii - 1);
							ref_ = ref1_.substr(ii-1, 1);
							alt_ = alt1_.substr(ii-1, jj-ii+1);
						} else {
							varPos = std::to_string(stoul(pos) + ii);
							ref_ = ref1_.substr(ii, 1);
							alt_ = alt1_.substr(ii, jj-ii);
						}
						alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor1);
						posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
						inpos += (jj-ii);
						ii = jj + 1;
					} else if ( ref1_[ii] != alt1_[ii] && ii+1 < max && alt1_.substr(ii+1,1) == string{"-"} ) {  // multiple sites deletion variation
						jj = ii + 1;
						while ( jj+1 < max && alt1_.substr(jj+1,1) == string{"-"} ) jj++;
						if ( alt1_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii - 1);
							ref_ = ref1_.substr(ii-1, jj-ii+1);
							alt_ = alt1_.substr(ii-1, 1);
						} else {
							varPos = std::to_string(stoul(pos) + ii);
							ref_ = ref1_.substr(ii, jj-ii);
							alt_ = alt1_.substr(ii, 1);
						}
						alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor1);
					   posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
						ii = jj + 1;
					} else if ( ref1_[ii] != alt1_[ii] ) {
						if ( ref1_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii - 1);
							ref_ = ref1_.substr(ii-1, 1);
							alt_ = alt1_.substr(ii-1, 2);
						} else if ( alt1_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii -1);
							ref_ = ref1_.substr(ii-1, 2);
							alt_ = alt1_.substr(ii-1, 1);
						} else {
							varPos = std::to_string(stoul(pos) + ii);
							ref_ = ref1_.substr(ii, 1);
							alt_ = alt1_.substr(ii, 1);
						}
						varPos = std::to_string(stoul(pos) + ii);
						alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor1);
						posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
						ii++;
					} else {
						ii++;
					}
				}
				// alignment the other alleles
				max = ref2_.size() >= alt2_.size() ? ref2_.size() : alt2_.size();
				min = ref2_.size() <= alt2_.size() ? ref2_.size() : alt2_.size();
				for ( auto j = max-1; j >= 0; j-- ) {
					if ( ref2_.size() == alt2_.size() ) {
						break;
					} else if ( ref2_.size() < alt2_.size() ) {
						auto coor = min - ( max - j );
						if ( ref2_[coor] != alt2_[j] ) {
							ref2_.insert(coor+1, "-");
							min++;
						}
					} else {
						auto coor = min - ( max - j );
						if ( alt2_[coor] != ref2_[j] ) {
							alt2_.insert(coor+1, "-");
							min++;
						}
					}
				}
				// std::cout << ref1_ << "\n" << alt1_ << "\n" << ref2_ << "\n" << alt2_ << std::endl;
				ii = 0; jj = 1; inpos = 0;
				while ( ii < max ) {
					if ( ref2_[ii] != alt2_[ii] && ii+1 < max && ref2_.substr(ii+1,1) == string{"-"} ) {         // multiple sites insert variation
						jj = ii + 1;
						while ( jj+1 < max && ref2_.substr(jj+1,1) == string{"-"} ) jj++;
						if ( ref2_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii - 1);
							ref_ = ref2_.substr(ii-1, 1);
							alt_ = alt2_.substr(ii-1, jj-ii+1);
						} else {
							varPos = std::to_string(stoul(pos) + ii);
							ref_ = ref2_.substr(ii, 1);
							alt_ = alt2_.substr(ii, jj-ii);
						}
						if ( posVarMap.find(varPos) == posVarMap.end() ) {
								alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor2);
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );

						} else {
							if ( posVarMap.at(varPos)[0] == make_tuple(varPos, ref_, alt_, infor1) ) {
								alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor);
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
							} else {
								alterAllelesVector = posVarMap.at(varPos);
								alterAllelesVector.push_back( make_tuple(varPos, ref_, alt_, infor2) );
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
								alterAllelesVector.erase(alterAllelesVector.begin()+1);
							}
						}
						inpos += (jj-ii);
						ii = jj + 1;
					} else if ( ref2_[ii] != alt2_[ii] && ii+1 < max && alt2_.substr(ii+1,1) == string{"-"} ) {  // multiple sites deletion variation
						jj = ii + 1;
						while ( jj+1 < max && alt2_.substr(jj+1,1) == string{"-"} ) jj++;
						if ( alt2_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii - 1);
							ref_ = ref2_.substr(ii-1, jj-ii+1);
							alt_ = alt2_.substr(ii-1, 1);
						} else {
							varPos = std::to_string(stoul(pos) + ii);
							ref_ = ref2_.substr(ii, jj-ii);
							alt_ = alt2_.substr(ii, 1);
						}
						if ( posVarMap.find(varPos) == posVarMap.end() ) {
							alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor2);
							posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
						} else {
							if ( posVarMap.at(varPos)[0] == make_tuple(varPos, ref_, alt_, infor1) ) {
								alterAllelesVector[0] = make_tuple(varPos, ref2_.substr(ii, jj-ii), alt2_.substr(ii, 1), infor);
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
							} else {
								alterAllelesVector = posVarMap.at(varPos);
								alterAllelesVector.push_back( make_tuple(varPos, ref_, alt_, infor2) );
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
								alterAllelesVector.erase(alterAllelesVector.begin()+1);
							}
						}
						ii = jj + 1;
					} else if ( ref2_[ii] != alt2_[ii] ) {
						if ( ref2_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii - 1);
							ref_ = ref2_.substr(ii-1, 1);
							alt_ = alt2_.substr(ii-1, 2);
						} else if ( alt2_.substr(ii, 1) == "-" ) {
							varPos = std::to_string(stoul(pos) + ii -1);
							ref_ = ref2_.substr(ii-1, 2);
							alt_ = alt2_.substr(ii-1, 1);
						} else {
							varPos = std::to_string(stoul(pos) + ii);
							ref_ = ref2_.substr(ii, 1);
							alt_ = alt2_.substr(ii, 1);
						}
						if ( posVarMap.find(varPos) == posVarMap.end() ) {
							alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor2);
							posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
						} else {
							if ( posVarMap.at(varPos)[0] == make_tuple(varPos, ref_, alt_, infor1) ) {
								alterAllelesVector[0] = make_tuple(varPos, ref_, alt_, infor);
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
							} else {
								alterAllelesVector = posVarMap.at(varPos);
								alterAllelesVector.push_back( make_tuple(varPos, ref_, alt_, infor2) );
								posVarMap.insert( std::make_pair(varPos, alterAllelesVector) );
								alterAllelesVector.erase(alterAllelesVector.begin()+1);
							}
						}
						ii++;
					} else {
						ii++;
					}
				}
				// travel all variation elements
				for ( auto const & posKey : posVarMap ) {
					for ( auto const & var : posKey.second ) {
						std::cout << std::get<0>(var) << " " << std::get<1>(var) << " " << std::get<2>(var) << std::endl;
						knownVar = getRelatedRecords(chr,std::get<0>(var),std::get<1>(var),std::get<2>(var),recordDb, false);
						if ( knownVar.size() >= 1 ) {
							writtenAllRecord(anno_file, std::get<3>(var), freqDb, knownVar);
							std::cout << "One end" << std::endl;
						} else {      // denovo variation site
							string novoSite = chr+"\t"+std::get<0>(var)+"\t"+std::get<1>(var)+"\t"+std::get<2>(var);
							DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
							if ( 0 == unknownSites.invoke_command(command) ) {
								knownVar = unknownSites.readResult();
								writtenAllRecord(anno_file, std::get<3>(var), freqDb, knownVar);
								insertToDatabase(mysqlConn, knownVar, recordDb.getOtherAnnoFields());
							} else {
								Logger.error("Unannotated denovo site: " + chr + " " + pos + " " + ref + " " + alt);
							}
						}
					}
				}
			}
		}
	}
}
