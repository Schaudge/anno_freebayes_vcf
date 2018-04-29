/* 
 *  Make a de novo small variation record analysis.
 *  Written By Schaudge King
 *  Date: 2017-01-20
 */
#ifndef DE_INVOKE_H
#define DE_INVOKE_H
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>

class DenovoAnalysis {

public:
	DenovoAnalysis( const std::string & inDir_, const std::string & annoInput_, const std::string & rawVar_ ) 
: inDir(inDir_), annoInput(annoInput_), rawVar(rawVar_)
{
		auto i = rawVar.find_last_of("\t");
		std::string alt = rawVar.substr(i+1), prefix = rawVar.substr(0,i);
		i = prefix.find_last_of("\t");
		std::string ref = prefix.substr(i+1);
		prefix = prefix.substr(0,i);
		i = prefix.find_last_of("\t");
		std::string position = prefix.substr(i+1);
		prefix = prefix.substr(0,i+1);

		std::ofstream annoFile( inDir + "/" + annoInput );
		auto _ref = ref.size(), _alt = alt.size();
		if ( _ref == 1 && _alt == 1 ) {
			annoFile << prefix << position << "\t" << position << "\t" << ref << "\t" << alt << "\n";
		} else if ( _ref > 2 && ref[0] != alt[0] && ref[_ref-1] != alt[_alt-1] && ref.substr(1,_ref-2) == alt.substr(1,_alt-2) ) { 
			annoFile << prefix << position << "\t" << position << "\t" << ref[0] << "\t" << alt[0] << "\n";
			annoFile << prefix << stoul(position)+_ref-1 << "\t" << stoul(position)+_ref-1 << "\t" << ref[_ref-1] << "\t" << alt[_alt-1] << "\n";
		} else {
			unsigned long j = 0;
			unsigned long max = _ref > _alt ? _alt : _ref;
			for ( j=0; j<max; ++j ) {
				if ( ref[ref.size()-1-j] != alt[alt.size()-1-j] ) {
					ref.erase(_ref-j); alt.erase(_alt-j);
					break;
				}
			}
			//std::cout << "j= "<< j << " "<< ref << " " << alt << "\n";
			if ( ref.substr(0,max-j) == alt.substr(0,max-j) ) {
				ref.erase(0,max-j); alt.erase(0,max-j);
				if ( ref == "" ) ref = "-"; else if ( alt == "" ) alt = "-";
			} else {
				for ( i=0; i<=(max-j); ++i ) {
					if ( ref[i] != alt[i] ) {
						ref.erase(0,i); alt.erase(0,i);
						break;
					}
				}
			}
			//std::cout << i << " " << ref << " " << alt << "\n";
			annoFile << prefix << stoul(position)+i << "\t" << stoul(position)+i+ref.size()-1 << "\t" << ref << "\t" << alt << "\n";
		}
};
	int invoke_command( const std::string & annoBinComm ) {
		std::string command_line = annoBinComm + " " + inDir + "/" + annoInput;
		if ( 0 == system( command_line.c_str() ) ) {
			return 0;
		} else {
			return -1;
		}
	};
	std::vector< std::tuple<std::string,std::string> > readResult() {
		std::ifstream resultFile( "./denovo.hg38_multianno.new" );
		std::string annoRecord;
		std::vector< std::tuple<std::string,std::string> >  varAnnoResult;
		while ( getline(resultFile, annoRecord) ) {
			auto split = annoRecord.find_first_of("|");
			if ( split != std::string::npos ) {
				varAnnoResult.push_back( make_tuple(annoRecord.substr(0,split), annoRecord.substr(split+1)) );
			}
		}
		return varAnnoResult;
	}
private:
	const std::string inDir;
	const std::string annoInput;
	const std::string rawVar;
};
#endif

