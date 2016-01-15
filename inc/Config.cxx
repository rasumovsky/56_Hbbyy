////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  A brief class that handles reading input from a configuration text file   //
//  based on root's TEnv format, i.e. key-value pairs, and shamelessly stolen //
//  from the HGamAnalysisTools package.                                       //
//                                                                            //
//  Author: Dag Gillberg                                                      //
//  Appropriator: Andrew Hard                                                 //
//  Email: ahard@cern.ch                                                      //
//  Date: 03/08/2015                                                          //
//                                                                            //
//  Config:                                                                   //
//    class to read settings from text files                                  //
//    relies on root's TEnv                                                   //
//                                                                            //
//  Usage:                                                                    //
//    Config settings("Hgamma.config");                                       //
//    TString gamContainerName = settings.getStr("PhotonContainer");          //
//    TString elContainerName  = settings.getStr("ElectronContainer");        //
//    vector<TString> systShifts = settings.getStrV("Systematics");           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "Config.h"

/**
   Class constructor.
*/
Config::Config(TString fileName) {
  addFile(fileName);
}

/**
   Ensures that there is a value in the database assocated with key if not, 
   abort with error message
*/
void Config::ensureDefined(TString key) {
  if (!isDefined(key)) {
    std::cout << "Config: no value found for " << key << std::endl;
    exit(0);
  }
}

/**
   returns true if the key is defined
*/
bool Config::isDefined(TString key) {
  return m_env.Defined(key);
}

/**
   Access a string value from the config database. Exception thrown if no entry exist.
*/
TString Config::getStr(TString key, bool expand) {
  ensureDefined(key);
  if (!expand) return m_env.GetValue(key,"");
  else return gSystem->ExpandPathName(m_env.GetValue(key,""));
}

/**
   Access a string value from the config database. Default value used if no entry exist.
*/
TString Config::getStr(TString key, TString dflt) {
  return m_env.GetValue(key, dflt);
}

/**
   Access an integer from the config database. Default value used if no entry exist.
*/
int Config::getInt(TString key) {
  ensureDefined(key);
  return m_env.GetValue(key,-99);
}

/**
   Access an integer from the config database. Exception thrown if no entry exist.
*/
int Config::getInt(TString key, int dflt) {
  return m_env.GetValue(key,dflt);
}

/**
   Access a boolean from the config database. Default value used if no entry exist
*/
bool Config::getBool(TString key, bool dflt) {
  return m_env.GetValue(key,dflt);
}

/**
   Access a boolean from the config database. Exception thrown if no entry exist.
*/
bool Config::getBool(TString key) {
  ensureDefined(key); return getBool(key,false);
}

/**
   Access a real number from the config database
*/
double Config::getNum(TString key) {
  ensureDefined(key);
  return m_env.GetValue(key,-99.0);
}

/**
   Access a vector of doubles from the config database
*/
double Config::getNum(TString key, double dflt) {
  return m_env.GetValue(key, dflt);
}

/**
   Access a vector of strings from the config database
*/
std::vector<TString> Config::getStrV(TString key) {
  ensureDefined(key);
  return vectorize(m_env.GetValue(key,"")," \t");
}

std::vector<double> Config::getNumV(TString key) {
  ensureDefined(key);
  return vectorizeNum(m_env.GetValue(key,"")," \t");
}

/**
   Prints the TEnv database to screen.
*/
void Config::printDB() {
  TIter next(m_env.GetTable());
  while (TEnvRec *er = (TEnvRec*) next()) {
    printf("  %-60s%s\n", Form("%s:", er->GetName()), er->GetValue());
  }
}

/**
   Add more user specified settings.
*/
void Config::addFile(TString fileName) {
  TString path(fileName);
  if (!fileExist(path) || path == "") {
    std::cout << "Cannot find settings file " << fileName
	      << "\n  also searched in " << path << std::endl;
    exit(0);
  }
  // settings read in by files should not overwrite values set by setValue()
  TEnv env;
  int status = env.ReadFile(path.Data(),EEnvLevel(0));
  if (status != 0) {
    std::cout << "Cannot read settings file " << fileName << std::endl;
    exit(0);
  }
  TIter next(env.GetTable());
  while (TEnvRec *er = (TEnvRec*) next()) {
    if (!isDefined(er->GetName())) {
      setValue(er->GetName(), er->GetValue());
    }
  }
}

/**
   Set value
*/
void Config::setValue(TString key, TString value) {
  m_env.SetValue(key,value);
}

std::vector<TString> Config::vectorize(TString str, TString sep) {
  std::vector<TString> result;
    TObjArray *strings = str.Tokenize(sep.Data());
    if (strings->GetEntries()==0) { delete strings; return result; }
    TIter istr(strings);
    while (TObjString* os=(TObjString*)istr()) {
      // the number sign and everything after is treated as a comment
      if (os->GetString()[0]=='#') break;
      result.push_back(os->GetString());
    }
    delete strings;
    return result;
  }
  
  // convert a text line containing a list of numbers to a vector<double>
std::vector<double> Config::vectorizeNum(TString str, TString sep) {
  std::vector<double> result; std::vector<TString> vecS = vectorize(str,sep);
    for (uint i=0;i<vecS.size();++i)
      result.push_back(atof(vecS[i]));
    return result;
  }

// checks if a given file or directory exist
bool Config::fileExist(TString fn) {
  return !(gSystem->AccessPathName(fn.Data()));
}
