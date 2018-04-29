/* 
 *  Make a snp records table for reference.
 *  Written By Schaudge King
 *  Date: 2017-01-20
 */
#ifndef REF_SNP_TAB_H
#define REF_SNP_TAB_H
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>

class RefSnpTable {

	friend std::vector<std::tuple<std::string, std::string> > getRelatedRecords( const std::string &, const std::string &,
			const std::string &, const std::string &,
			const RefSnpTable &, std::size_t & );

	friend std::vector<std::tuple<std::string, std::string> > getRelatedRecords( const std::string &, const std::string &,
			const std::string &, const std::string &,
			RefSnpTable &, bool flagMovable);
public:
	using RECORD_TYPE = std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, std::string>;

	RefSnpTable ( sql::Connection * & );
	RefSnpTable( std::ifstream & );
	void display_all();
	std::string getOtherAnnoFields() { return suffix_fields; };
	void display_region( const std::string & );
	const std::string getCurrentElement();

private:
	std::vector< std::shared_ptr<RECORD_TYPE> > refTable;
	std::string suffix_fields;

public:
	decltype(refTable)::const_iterator ref_flag;
};
#endif

std::vector<std::tuple<std::string, std::string> > getRelatedRecords( const std::string & chr, const std::string & pos, 
		const std::string & ref, const std::string & alt,
		const RefSnpTable & annoTable, std::size_t & ref_flag);
std::vector<std::tuple<std::string, std::string> > getRelatedRecords( const std::string &, const std::string &, 
		const std::string &, const std::string &,
		RefSnpTable &, bool flagMovable);
