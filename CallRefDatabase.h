/*
 * CallRefDatabase.h
 *
 *  Created on: 2018年4月27日
 *      Author: yuanshenran
 */

#ifndef SRC_ANNO_CALLREFDATABASE_H_
#define SRC_ANNO_CALLREFDATABASE_H_
#include <string>
#include <cppconn/driver.h>
#include <cppconn/statement.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>

class CallRefDatabase {

public:

	CallRefDatabase (const std::string & ipAddress, const std::string & user,
			const std::string & passwd, const std::string & database) :
		ipAddress_(ipAddress), user_(user), passwd_(passwd), database_(database) {
		mysql_driver = get_driver_instance();
		mysql_conn = mysql_driver->connect(ipAddress_, user_, passwd_);
		mysql_conn->setSchema(database_);
	}

	~ CallRefDatabase() {
		delete mysql_conn;
	}

	sql::Connection * & getDbConnection() {
		return mysql_conn;
	}


private:

	// created database connection
	sql::Driver * mysql_driver;
	sql::Connection * mysql_conn;

	// needed database configure
   std::string ipAddress_;
   std::string user_;
   std::string passwd_;
   std::string database_;

};


#endif /* SRC_ANNO_CALLREFDATABASE_H_ */
