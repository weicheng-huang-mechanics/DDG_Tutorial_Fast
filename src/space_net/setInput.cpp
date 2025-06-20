#include "setInput.h"

using namespace std;

setInput::setInput()
{
	AddOption("render", "visualization", render);
	AddOption("saveData", "should results be saved", saveData);
  AddOption("deltaTime", "deltaTime", deltaTime);
  AddOption("totalTime", "totalTime", totalTime);
  AddOption("YoungM", "YoungM", YoungM);
  AddOption("density", "density", density);
  AddOption("Possion", "Possion", Possion);
  AddOption("stol", "stol", stol);
  AddOption("forceTol", "forceTol", forceTol);
  AddOption("scaleRendering", "scaleRendering", scaleRendering);
  AddOption("maxIter", "maxIter", maxIter);
  AddOption("gVector", "gVector", gVector);
  AddOption("viscosity", "viscosity", viscosity);
  AddOption("deltaLength", "deltaLength", deltaLength);
  AddOption("rodRadius", "rodRadius", rodRadius);
  AddOption("boxSize", "boxSize", boxSize);
  AddOption("stiffness", "stiffness", stiffness);
  AddOption("dBar", "dBar", dBar);
  AddOption("dBarBody", "dBarBody", dBarBody);
  AddOption("stiffnessBody", "stiffnessBody", stiffnessBody);
}

setInput::~setInput()
{
	;
}

Option* setInput::GetOption(const string& name)
{
  if (m_options.find(name) == m_options.end()) 
  {
    cerr << "Option " << name << " does not exist" << endl;
  }
  return &(m_options.find(name)->second);
}

bool& setInput::GetBoolOpt(const string& name)
{
  return GetOption(name)->b;
}

int& setInput::GetIntOpt(const string& name)
{
  return GetOption(name)->i;
}

double& setInput::GetScalarOpt(const string& name)
{
  return GetOption(name)->r;
}

Vector3d& setInput::GetVecOpt(const string& name)
{
  return GetOption(name)->v;
}

string& setInput::GetStringOpt(const string& name)
{
  return GetOption(name)->s;
}

int setInput::LoadOptions(const char* filename)
{
  ifstream input(filename);
  if (!input.is_open()) 
  {
    cerr << "ERROR: File " << filename << " not found" << endl;
    return -1;
  }

  string line, option;
  istringstream sIn;
  string tmp;
  for (getline(input, line); !input.eof(); getline(input, line)) 
  {
    sIn.clear();
    option.clear();
    sIn.str(line);
    sIn >> option;
    if (option.size() == 0 || option.c_str()[0] == '#') continue;
    OptionMap::iterator itr;
    itr = m_options.find(option);
    if (itr == m_options.end()) 
    {
      cerr << "Invalid option: " << option << endl;
      continue;
    }
    if (itr->second.type == Option::BOOL) 
    {
      sIn >> tmp;
      if (tmp == "true" || tmp == "1") itr->second.b = true;
      else if (tmp == "false" || tmp == "0") itr->second.b = false;
    } 
    else if (itr->second.type == Option::INT) 
    {
      sIn >> itr->second.i;
    } 
    else if (itr->second.type == Option::DOUBLE) 
    {
      sIn >> itr->second.r;
    } 
    else if (itr->second.type == Option::VEC) 
    {
      Vector3d& v = itr->second.v;
      sIn >> v[0];
      sIn >> v[1];
      sIn >> v[2];
    } 
    else if (itr->second.type == Option::STRING) 
    {
      sIn >> itr->second.s;
    } else 
    {
      cerr << "Invalid option type" << endl;
    }
  }
  input.close();

  return 0;
}

int setInput::LoadOptions(int argc, char** argv)
{
  string option, tmp;
  int start = 0;
  while (start < argc && string(argv[start]) != "--") ++start;
  for (int i = start + 1; i < argc; ++i) 
  {
    option = argv[i];
    OptionMap::iterator itr;
    itr = m_options.find(option);
    if (itr == m_options.end()) 
    {
      cerr << "Invalid option on command line: " << option << endl;
      continue;
    }
    if (i == argc - 1) 
    {
      cerr << "Too few arguments on command line" << endl;
      break;
    }
    if (itr->second.type == Option::BOOL) 
    {
      tmp = argv[i+1]; ++i;
      if (tmp == "true" || tmp == "1") itr->second.b = true;
      if (tmp == "false" || tmp == "0") itr->second.b = false;
    } 
    else if (itr->second.type == Option::INT) 
    {
      itr->second.i = atoi(argv[i+1]); ++i;
    } 
    else if (itr->second.type == Option::DOUBLE) 
    {
      itr->second.r = atof(argv[i+1]); ++i;
    } 
    else if (itr->second.type == Option::VEC) 
    {
      if (i >= argc - 3) 
      {
        cerr << "Too few arguments on command line" << endl;
        break;
      }
      Vector3d& v = itr->second.v;
      v[0] = atof(argv[i+1]); ++i;
      v[1] = atof(argv[i+1]); ++i;
      v[2] = atof(argv[i+1]); ++i;
    } 
    else if (itr->second.type == Option::STRING) 
    {
      itr->second.s = argv[i+1]; ++i;
    } 
    else 
    {
      //cerr << "Invalid option type" << endl;
    }
  }
  return 0;
}
