#ifndef CONFIGREADER
#define CONFIGREADER

#include <string>

class ConfigReader{
  protected:
    std::string _filename;
  public:
    explicit ConfigReader(std::string filename);
    void setFilename(const std::string filename);
    std::string getFilename() const;
    virtual void setUp() = 0;
};

#endif
