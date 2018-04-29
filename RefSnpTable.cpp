/* 
 *  Make a snp records table for reference.
 *  Written By Schaudge King
 *  Date: 2017-01-20
 */
#include "RefSnpTable.h"
#include "DenovoInvoke.h"
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

RefSnpTable::RefSnpTable ( sql::Connection * & mysqlConn)
{
	sql::Statement * stmt = mysqlConn->createStatement();
	sql::ResultSet * res = stmt->executeQuery("SELECT * FROM nanfang_gene ORDER BY chr,pos ASC");
	sql::ResultSetMetaData * res_meta = res->getMetaData();
	int colCounts = res_meta->getColumnCount();
	for (int i = 8; i <= colCounts; i++) {
		suffix_fields += ( "\t" + res_meta->getColumnLabel(i) );
	}
	string line, chr, pos, ref, alt, rid, sid, others;
	while (res->next())
	{
		rid = res->getString("search_id");
		sid = res->getString("snp_id");
		chr = res->getString("chr");
		pos = res->getString("pos");
		ref = res->getString("ref");
		alt = res->getString("alt");
		others = "";
		for (int i = 8; i <= colCounts; i++) {
			others += ( "\t" + res->getString(i) );
		}
		shared_ptr< RECORD_TYPE > record_ptr = make_shared <RECORD_TYPE> ( make_tuple(chr, pos, ref, alt, rid, sid, others) );
		refTable.push_back(record_ptr);
	}
	ref_flag = refTable.cbegin();
	delete stmt;
	delete res;
}

RefSnpTable::RefSnpTable (ifstream & record_file)
{
	suffix_fields = "\ttype\tgene_symbol\tfunctional_consequence_dbsnp147\tclinical_significance_dbsnp147\taa_change_dbsnp147\tglobal_MAF_dbsnp147\tClinicalSignificance_ClinVar\tother_id_ClinVar\tMAF_1000g\tEAS_1000g\tprimary_site_cosmic\tFATHMM_prediction_cosmic\tsift_score_ljb26\tpolyphen2_HDIV_score_ljb26\tfreq_esp6500\ttag_hgmd\tdisease_hgmd\tpmid_hgmd\tClinically_Importance_BIC\tCategory_BIC\tdbSNP147_china_freq\tsift_score_ensembl\tsift_class_ensembl\tpolyphen_score_ensembl\tpolyphen_class_ensembl";
	string line, chr, pos, ref, alt, rid, sid, others;
	while ( record_file && getline(record_file, line) ) {
		istringstream record(line);
		record >> rid; record >> sid;
		record >> chr; record >> pos;
		record >> ref; record >> alt;
		getline(record, others);
		shared_ptr< RECORD_TYPE > record_ptr = make_shared <RECORD_TYPE> ( make_tuple(chr, pos, ref, alt, rid, sid, others) );
		refTable.push_back(record_ptr);
	}
	ref_flag = refTable.cbegin();
}

vector< tuple<string, string> > getRelatedRecords( const string & chr, const string & pos, const string & ref, const string & alt, 
		const RefSnpTable & annoTable, size_t & ref_flag)
{
	auto iter = annoTable.refTable.cbegin() + ref_flag, iter_end = annoTable.refTable.cend();
	auto position = stoul(pos), end = position + ref.size();
	vector< tuple<string,string> > identicalRecords;
	vector< decltype(iter) > relatedSet;
	string ref_ = ref;
	for ( ; iter != iter_end; ++iter) {
		++ref_flag;
		ref_ = ref;
		string rchr = get<0>(**iter);
		if ( rchr.compare(chr) == 0 ) {
			auto rpos = stoul(get<1>(**iter));
			if ( rpos >= position ) {
				auto rref = get<2>(**iter);
				if ( rpos + rref.size() <= end ) {
					auto ralt = get<3>(**iter);
					if ( ref.size() == 1 && alt.size() == 1 ) {   // case snp
						if ( rref != "-" && ralt == alt ) {
							identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
							return identicalRecords;
						}
					} else if ( ref.size() < alt.size() ) {       // case insert, possibality indel
						if ( rref == "-" ) {
							if ( ref_.insert(rpos-position+1, ralt) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							} else if ( ralt != "-" ) {
								relatedSet.push_back( iter );
							}
						} else if ( rref.size() > 1 || ralt.size() > 1 ) {      // case indel
							if ( ref_.replace(rpos-position, rref.size(), ralt) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							}
						}
					} else {                                      // case delete, possibality indel
						if ( ralt == "-" ) {
							if ( ref_.erase(rpos-position, rref.size() ) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							} else if ( rref != "-" ) {
								relatedSet.push_back( iter );
							}
						} else if ( rref.size() > 1 || ralt.size() > 1 ) {      // case indel
							if ( ref_.replace(rpos-position, rref.size(), ralt) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							}
						}
					}
				} else if ( rpos >= end ) {
					break;
				}
			}
		} else if ( rchr.compare(chr) > 0 ) {
			break;
		}
	}
	/* Travel all two elements' combination for complex insert/delete, from larger position to lower position.
	          rpos2 <----- rpos1  */
	if ( relatedSet.size() >= 2 ) {
		string ref1_, ref2_;
		auto combin2 = relatedSet.size() - 1, accum2 = relatedSet.size() - 2;
		size_t rpos1=0, rpos2=0;
		string rref1, ralt1, rref2, ralt2;
		while ( combin2 >= 1 ) {
			ref1_ = ref;
			rpos1 = stoul(get<1>( **(relatedSet[combin2]) ));
			rref1 = get<2>( **(relatedSet[combin2]) );
			ralt1 = get<3>( **(relatedSet[combin2]) );
			if ( rref1 == "-" ) {
				ref1_.insert( rpos1-position+1, ralt1 );
			} else if ( ralt1 == "-" ) {
				ref1_.erase( rpos1-position, rref1.size() );
			} else {
				ref1_.replace( rpos1-position, rref1.size(), ralt1 );
			}
			for ( accum2=combin2-1; accum2 > 0 ; --accum2 ) {
				ref2_ = ref1_;
				rpos2 = stoul(get<1>( **(relatedSet[accum2]) ));
				rref2 = get<2>( **(relatedSet[accum2]) );
				ralt2 = get<3>( **(relatedSet[accum2]) );
				if ( rref2 == "-" ) {
					ref2_.insert( rpos2-position+1, ralt2 );
				} else if ( ralt2 == "-" ) {
					ref2_.erase( rpos2-position, rref2.size() );
				} else {
					ref2_.replace( rpos2-position, rref2.size(), ralt2 );
				}
				if ( 0 == alt.compare(ref2_) ) {
					identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[combin2]))+"\t"+get<5>(**(relatedSet[combin2]))+"\t"+chr+"\t"
							+ get<1>(**(relatedSet[combin2]))+"\t"+rref1+"\t"+ralt1, get<6>(**(relatedSet[combin2]))) );
					identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[accum2]))+"\t"+get<5>(**(relatedSet[accum2]))+"\t"+chr+"\t"
							+ get<1>(**(relatedSet[accum2]))+"\t"+rref2+"\t"+ralt2, get<6>(**(relatedSet[accum2]))) );
					return identicalRecords;
				}
			}
			rpos2 = stoul(get<1>( **(relatedSet[0]) ));
			rref2 = get<2>( **(relatedSet[0]) );
			ralt2 = get<3>( **(relatedSet[0]) );
			if ( rref2 == "-" ) {
				ref1_.insert( rpos2-position+1, ralt2 );
			} else if ( ralt2 == "-" ) {
				ref1_.erase( rpos2-position, rref2.size() );
			} else {
				ref1_.replace( rpos2-position, rref2.size(), ralt2 );
			}
			if ( 0 == alt.compare(ref1_) ) {
				identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[combin2]))+"\t"+get<5>(**(relatedSet[combin2]))+"\t"+chr+"\t"
						+ get<1>(**(relatedSet[combin2]))+"\t"+rref1+"\t"+ralt1, get<6>(**(relatedSet[combin2]))) );
				identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[accum2]))+"\t"+get<5>(**(relatedSet[accum2]))+"\t"+chr+"\t"
						+ get<1>(**(relatedSet[accum2]))+"\t"+rref2+"\t"+ralt2, get<6>(**(relatedSet[accum2]))) );
				return identicalRecords;
			}
			--combin2;
		}
	}
	return identicalRecords;
}


void RefSnpTable::display_all() {
	auto iter = refTable.cbegin(), iter_end = refTable.cend();
	string rchr, rpos, rref, ralt, snpId, other;
	auto count = 0;
	for ( ; iter != iter_end; ++iter) {
		rchr = get<0>( **iter );
		rpos = get<1>( **iter );
		rref = get<2>( **iter );
		ralt = get<3>( **iter );
		snpId = get<4>( **iter );
		other = get<6>( **iter );
		cout << "Record: " << snpId << " Infor: " << rchr << " " << rpos << " " << rref << " " << ralt << " " << other << "\n";
		count++;
	}
	cout << "----------------------------------------------\n";
	cout << "Total records counts are " << count << "." << endl;
}

const string RefSnpTable::getCurrentElement() {
	auto iter = RefSnpTable::ref_flag;
	string rchr, rpos, rref, ralt, snpId, other;
	rchr = get<0>( **iter );
	rpos = get<1>( **iter );
	rref = get<2>( **iter );
	ralt = get<3>( **iter );
	snpId = get<4>( **iter );
	other = get<6>( **iter );
	return rchr+" "+rpos+" "+rref+" "+ralt;
}

void RefSnpTable::display_region( const string & region ) {
	auto i = region.find(":", 0);
	string chr = region.substr(0,i), coordinate = region.substr(i+1);
	i = coordinate.find("-");
	std::size_t start = stoi(coordinate.substr(0,i)), end = stoi(coordinate.substr(i+1));
	cout << start << " " << end << endl;
	auto iter = refTable.cbegin(), iter_end = refTable.cend();
	string rchr, rpos, rref, ralt, snpId, other;
	std::size_t count = 0;
	for ( ; iter != iter_end; ++iter) {
		rchr = get<0>(**iter);
		if ( rchr.compare(chr) == 0 ) {
			std::size_t rpos = stoi(get<1>(**iter));
			if ( rpos >= start ) {
				auto rref = get<2>(**iter);
				if ( rpos + rref.size() <= end ) {
					ralt = get<3>( **iter );
					snpId = get<4>( **iter );
					other = get<6>( **iter );
					cout << "Record: " << snpId << " Infor: " << rchr << " " << rpos << " " << rref << " " << ralt << " " << other << "\n";
					count++;
				} else if ( rpos > end ) {
					break;
				}
			}
		} else if ( rchr.compare(chr) > 0 ) {
			break;
		}
	}
	cout << "----------------------------------------------\n";
	cout << "Total records counts are " << count << "." << endl;
}


vector< tuple<string, string> > getRelatedRecords( const string & chr, const string & pos, const string & ref, const string & alt, 
		RefSnpTable & annoTable, bool flagMovable=true)
{
	auto iter = annoTable.ref_flag, iter_end = annoTable.refTable.cend();
	auto position = stoul(pos), end = position + ref.size();
	vector< tuple<string,string> > identicalRecords;
	vector< decltype(iter) > relatedSet;
	string ref_ = ref;
	for ( ; iter != iter_end; ++iter) {
		if (flagMovable) ++(annoTable.ref_flag);
		ref_ = ref;
		string rchr = get<0>(**iter);
		if ( rchr.compare(chr) == 0 ) {
			auto rpos = stoul(get<1>(**iter));
			if ( rpos >= position ) {
				auto rref = get<2>(**iter);
				if ( rpos + rref.size() <= end ) {
					auto ralt = get<3>(**iter);
					if ( ref.size() == 1 && alt.size() == 1 ) {   // case snp
						if ( rref != "-" && ralt == alt ) {
							identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
							return identicalRecords;
						}
					} else if ( ref.size() < alt.size() ) {       // case insert, possibality indel
						if ( rref == "-" ) {
							if ( ref_.insert(rpos-position+1, ralt) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							} else if ( ralt != "-" ) {
								relatedSet.push_back( iter );
							}
						} else if ( rref.size() > 1 || ralt.size() > 1 ) {      // case indel
							if ( ref_.replace(rpos-position, rref.size(), ralt) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							}
						}
					} else {                                      // case delete, possibality indel
						if ( ralt == "-" ) {
							if ( ref_.erase(rpos-position, rref.size() ) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							} else if ( rref != "-" ) {
								relatedSet.push_back( iter );
							}
						} else if ( rref.size() > 1 || ralt.size() > 1 ) {      // case indel
							if ( ref_.replace(rpos-position, rref.size(), ralt) == alt ) {
								identicalRecords.push_back( make_tuple(get<4>(**iter)+"\t"+get<5>(**iter)+"\t"+rchr+"\t"+get<1>(**iter)+"\t"+rref+"\t"+ralt, get<6>(**iter)) );
								return identicalRecords;
							}
						}
					}
				} else if ( rpos >= end ) {
					break;
				}
			}
		} else if ( rchr.compare(chr) > 0 ) {
			break;
		}
	}
	/* Travel all two elements' combination for complex insert/delete, from larger position to lower position.
	          rpos2 <----- rpos1  */
	if ( relatedSet.size() >= 2 ) {
		string ref1_, ref2_;
		auto combin2 = relatedSet.size() - 1, accum2 = relatedSet.size() - 2;
		size_t rpos1=0, rpos2=0;
		string rref1, ralt1, rref2, ralt2;
		while ( combin2 >= 1 ) {
			ref1_ = ref;
			rpos1 = stoul(get<1>( **(relatedSet[combin2]) ));
			rref1 = get<2>( **(relatedSet[combin2]) );
			ralt1 = get<3>( **(relatedSet[combin2]) );
			if ( rref1 == "-" ) {
				ref1_.insert( rpos1-position+1, ralt1 );
			} else if ( ralt1 == "-" ) {
				ref1_.erase( rpos1-position, rref1.size() );
			} else {
				ref1_.replace( rpos1-position, rref1.size(), ralt1 );
			}
			for ( accum2=combin2-1; accum2 > 0 ; --accum2 ) {
				ref2_ = ref1_;
				rpos2 = stoul(get<1>( **(relatedSet[accum2]) ));
				rref2 = get<2>( **(relatedSet[accum2]) );
				ralt2 = get<3>( **(relatedSet[accum2]) );
				if ( rref2 == "-" ) {
					ref2_.insert( rpos2-position+1, ralt2 );
				} else if ( ralt2 == "-" ) {
					ref2_.erase( rpos2-position, rref2.size() );
				} else {
					ref2_.replace( rpos2-position, rref2.size(), ralt2 );
				}
				if ( 0 == alt.compare(ref2_) ) {
					identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[combin2]))+"\t"+get<5>(**(relatedSet[combin2]))+"\t"+chr+"\t"
							+ get<1>(**(relatedSet[combin2]))+"\t"+rref1+"\t"+ralt1, get<6>(**(relatedSet[combin2]))) );
					identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[accum2]))+"\t"+get<5>(**(relatedSet[accum2]))+"\t"+chr+"\t"
							+ get<1>(**(relatedSet[accum2]))+"\t"+rref2+"\t"+ralt2, get<6>(**(relatedSet[accum2]))) );
					return identicalRecords;
				}
			}
			rpos2 = stoul(get<1>( **(relatedSet[0]) ));
			rref2 = get<2>( **(relatedSet[0]) );
			ralt2 = get<3>( **(relatedSet[0]) );
			if ( rref2 == "-" ) {
				ref1_.insert( rpos2-position+1, ralt2 );
			} else if ( ralt2 == "-" ) {
				ref1_.erase( rpos2-position, rref2.size() );
			} else {
				ref1_.replace( rpos2-position, rref2.size(), ralt2 );
			}
			if ( 0 == alt.compare(ref1_) ) {
				identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[combin2]))+"\t"+get<5>(**(relatedSet[combin2]))+"\t"+chr+"\t"
						+ get<1>(**(relatedSet[combin2]))+"\t"+rref1+"\t"+ralt1, get<6>(**(relatedSet[combin2]))) );
				identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[accum2]))+"\t"+get<5>(**(relatedSet[accum2]))+"\t"+chr+"\t"
						+ get<1>(**(relatedSet[accum2]))+"\t"+rref2+"\t"+ralt2, get<6>(**(relatedSet[accum2]))) );
				return identicalRecords;
			}
			--combin2;
		}
	}
	/* If there is no identical combination of dbsnp database records to generate the variation result,
	 * we give a left alignment annotation result for normalization, see also bcftools !
	 * See https://genome.sph.umich.edu/wiki/Variant_Normalization#Algorithm_for_Normalization
	 * Reference: 1. https://academic.oup.com/bioinformatics/article/31/13/2202/196142
	 *            2. https://academic.oup.com/bioinformatics/article/33/7/964/2623048
	 */
	if ( flagMovable && relatedSet.size() > 0 ) {
		if (identicalRecords.size() < 1) {                 // Make ref and alt sitewise aligned
			string ref_ = ref, alt_ = alt;
			unsigned long max = ref.size() >= alt.size() ? ref.size() : alt.size();
			unsigned long min = ref.size() <= alt.size() ? ref.size() : alt.size();
			for ( auto j = max-1; j >= 0; j-- ) {
				if ( ref_.size() == alt_.size() ) {
					break;
				} else if ( ref_.size() < alt_.size() ) {
					auto coor = min - ( max - j );
					if ( ref_[coor] != alt_[j] ) {
						ref_.insert(coor+1, "-");
						min++;
					}
				} else {
					auto coor = min - ( max - j );
					if ( alt_[coor] != ref_[j] ) {
						alt_.insert(coor+1, "-");
						min++;
					}
				}
			}
			// search identical result in the relatedSet !!!
			unsigned long ii = 0, jj = 1, inpos = 0;
			while ( ii < max ) {
				bool findflag = false;
				if ( ref_[ii] != alt_[ii] && ii+1 < max && ref_.substr(ii+1,1) == string{"-"} ) {         // multiple sites insert variation
					jj = ii + 1;
					while ( jj+1 < max && ref_.substr(jj+1,1) == string{"-"} ) jj++;
					for (unsigned long k=0; k < relatedSet.size()-1; k++) {
						unsigned long rpos = stoul(get<1>(**(relatedSet[k])));
						string rref = get<2>(**(relatedSet[k]));
						string ralt = get<3>(**(relatedSet[k]));
						if ( position+ii == rpos+inpos && ref_.substr(ii,jj-ii) == rref && alt_.substr(ii,jj-ii) == ralt ) {
							identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[k]))+"\t"+get<5>(**(relatedSet[k]))+"\t"+chr+"\t"
									+ get<1>(**(relatedSet[k]))+"\t"+rref+"\t"+ralt, get<6>(**(relatedSet[k]))) );
							findflag = true;
							break;
						}
					}
					if ( ! findflag ) {
						string novoSite = chr+"\t"+to_string(position+ii)+"\t"+ref_.substr(ii,1)+"\t"+alt_.substr(ii,1);
						DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
						if ( 0 == unknownSites.invoke_command("/gcbi/analyze/integration/denovo.sh") ) {
							identicalRecords.push_back( unknownSites.readResult()[0] );
						}
					}
					inpos += (jj-ii);
					ii = jj + 1;
				} else if ( ref_[ii] != alt_[ii] && ii+1 < max && alt_.substr(ii+1,1) == string{"-"} ) {  // multiple sites deletion variation
					jj = ii + 1;
					while ( jj+1 < max && alt_.substr(jj+1,1) == string{"-"} ) jj++;
					for (unsigned long k=0; k < relatedSet.size()-1; k++) {
						unsigned long rpos = stoul(get<1>(**(relatedSet[k])));
						string rref = get<2>(**(relatedSet[k]));
						string ralt = get<3>(**(relatedSet[k]));
						if ( position+ii == rpos && ref_.substr(ii,jj-ii) == rref && alt_.substr(ii,jj-ii) == ralt ) {
							identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[k]))+"\t"+get<5>(**(relatedSet[k]))+"\t"+chr+"\t"
									+ get<1>(**(relatedSet[k]))+"\t"+rref+"\t"+ralt, get<6>(**(relatedSet[k]))) );
							findflag = true;
							break;
						}
					}
					if ( ! findflag ) {
						string novoSite = chr+"\t"+to_string(position+ii)+"\t"+ref_.substr(ii,1)+"\t"+alt_.substr(ii,1);
						DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
						if ( 0 == unknownSites.invoke_command("/gcbi/analyze/integration/denovo.sh") ) {
							identicalRecords.push_back( unknownSites.readResult()[0] );
						}
					}
					ii = jj + 1;
				} else if ( ref_[ii] != alt_[ii] ) {                                                      // single site variation
					for (unsigned long k=0; k < relatedSet.size()-1; k++) {
						unsigned long rpos = stoul(get<1>(**(relatedSet[k])));
						string rref = get<2>(**(relatedSet[k]));
						string ralt = get<3>(**(relatedSet[k]));
						if ( position+ii == rpos && ref_.substr(ii,1) == rref && alt_.substr(ii,1) == ralt ) {
							identicalRecords.push_back( make_tuple(get<4>(**(relatedSet[k]))+"\t"+get<5>(**(relatedSet[k]))+"\t"+chr+"\t"
									+ get<1>(**(relatedSet[k]))+"\t"+rref+"\t"+ralt, get<6>(**(relatedSet[k]))) );
							findflag = true;
							break;
						}
					}
					if ( ! findflag ) {
						string novoSite = chr+"\t"+to_string(position+ii)+"\t"+ref_.substr(ii,1)+"\t"+alt_.substr(ii,1);
						DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
						if ( 0 == unknownSites.invoke_command("/gcbi/analyze/integration/denovo.sh") ) {
							identicalRecords.push_back( unknownSites.readResult()[0] );
						}
					}
					ii++;
				}
				ii++;
			}
		}
	}
	return identicalRecords;
}

