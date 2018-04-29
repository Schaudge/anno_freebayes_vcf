/* 
 *  An unit test for class RefSnpTable.
 *  Written By Schaudge King
 *  Date: 2017-01-20
 */
#include "RefSnpTable.h"
#include "DenovoInvoke.h"
#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <fstream>
#include <iostream>

using namespace std;
								   
int main1( int argc, char * argv[] ) {

	vector<tuple<string, string> > result;
        /*
	DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", "17\t43097346\tTAAAAAAATGT\tTAAAAAATGT");
        auto iter = result.cbegin(), iter_end = result.cend();
        if ( 0 == unknownSites.invoke_command("/gcbi/analyze/integration/denovo.sh") ) {
            result = unknownSites.readResult();
            iter = result.cbegin(); iter_end = result.cend();
            for ( ; iter != iter_end; ++iter) {
                 cout << get<0>(*iter) << "\n";
            }
        } else {
            cout << "log this record...\n";
        }
        */
        cout << "------------------------------\n";

	//ifstream ref_snp_file("/home/wangxiaozhi/BioRef/chgc/BRCA_CARE_SNP147.index");
	//RefSnpTable recordTest = RefSnpTable(ref_snp_file);
	sql::Driver *mysql_driver;
	sql::Connection * mysql_con;
	/* Create a connection */
	mysql_driver = get_driver_instance();
	mysql_con = mysql_driver->connect("tcp://192.168.2.10:3306", "devuser", "111111");
	/* Connect to the MySQL database */
	mysql_con->setSchema("nanfang_gene");
	RefSnpTable recordTest = RefSnpTable(mysql_con);
	//recordTest.display_all();
        //recordTest.display_region("17:43093671-43093681");
	string novoSite = "13\t32333165\tTGGCC\tTGCA";
	DenovoAnalysis unknownSites = DenovoAnalysis("/gcbi/analyze/integration", "denovo_var.avinput", novoSite);
	vector<tuple<string,string> > knownVar = unknownSites.readResult();
	auto iter = knownVar.cbegin(), iter_end = knownVar.cend();
	for ( ; iter != iter_end; ++iter ) {
	     cout << get<0>(*iter) << get<1>(*iter) << "\n";			             
	}
	/*
	result = getRelatedRecords( "13", "32337326" , "A" , "G" , recordTest);
        iter = result.cbegin();
        iter_end = result.cend();
        for ( ; iter != iter_end; ++iter) {
             cout << get<0>(*iter) << "\n";
        }
	cout << recordTest.getCurrentElement() << endl;
        size_t start = 0;
        result = getRelatedRecords( "17", "43097346" , "TAAAAAAATGT" , "TAAAAAATGT" , recordTest, start);
        cout << "current reading records number: " << start << ".\n";
        iter = result.cbegin();
        iter_end = result.cend();
        for ( ; iter != iter_end; ++iter) {
             cout << get<0>(*iter) << "\n";
        }
        result = getRelatedRecords( "17", "43104058" , "AAAAAAAAAAAAGAAAAAAAAAAGAAAAGAAGAAGAAGAAGAAGAAGAAAACAAATGGT" , "AAAAAAAAAAAGAAAAGAAGAAGAAGAAGAAGAAAACAAATGGT" , recordTest, start);
        cout << "current reading records number: " << start << ".\n";
        iter = result.cbegin();
        iter_end = result.cend();
        for ( ; iter != iter_end; ++iter) {
             cout << get<0>(*iter) << "\n";
        }
	*/
	delete mysql_con;
	return 0;
}
