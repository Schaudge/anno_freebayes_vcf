/* 
 *  Make a frequence annotation table of China patient.
 *  Written By Schaudge King
 *  Date: 2017-03-15
 */
#ifndef FREQ_ANNO_TAB_H
#define FREQ_ANNO_TAB_H
#include <string>
#include <map>
#include <cppconn/statement.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>

class FreqAnnoTable {

public:

	FreqAnnoTable ( sql::Connection * & mysqlConn ) : column_names("\tartificial_grade\tfreq_total\tfreq_control\tfreq_case\tsignificance") {
		sql::Statement * stmt = mysqlConn->createStatement();
		sql::ResultSet * res = stmt->executeQuery("SELECT search_id,artificial_grade,freq_total,freq_control,freq_case,significance FROM variation_analyze_result");
		std::string rid, value;
		while (res->next())
		{
			rid = res->getString("search_id");
			value = "\t" + res->getString("artificial_grade") + "\t" + res->getString("freq_total")
					 + "\t" + res->getString("freq_control") + "\t" + res->getString("freq_case")
					 + "\t" + res->getString("significance");
			freqTable.emplace(rid, value);
		}
		delete stmt;
		delete res;
	};

	const std::string getOtherFields() const {
		return column_names;
	};

	std::map< std::string, std::string> * getFreqTable() {
		return &freqTable;
	};

	void display_all() {
		for (auto & pair : freqTable)
			std::cout << pair.first <<  ':' << pair.second << '\n';
	}

private:
	std::map< std::string, std::string> freqTable;
	std::string column_names;
};
#endif

