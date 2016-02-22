#ifndef CONFIGREADERANALYSIS
#define CONFIGREADERANALYSIS

#include "ConfigReader.h"

#include <string>

class ConfigReaderAnalysis : public ConfigReader{
  private:
    std::string _configName;
    bool _doEventLooper;
    long int _maxEvents;
		std::string _outputFile;
  public:
    ConfigReaderAnalysis(std::string filename);
    void setUp();
    void setConfigName(const std::string configname);
    void setDoEventLooper(const bool loop);
    void setMaxEvents(const long int i);
		void setOutputFile(std::string const outputfile);
    std::string getConfigName() const;
    bool getDoEventLooper() const;
    long int getMaxEvents() const;
};

#endif
