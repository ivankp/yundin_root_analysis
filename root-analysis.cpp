
#include <Python.h>
#include <cstdio>
#include <iostream>

#include "SelectorCommon.h"

#include "root-analysis.h"

using namespace RootAnalysis;

static char** g_argv = 0;

static std::vector<SelectorCommon*> selector_list;

bool RootAnalysis::Init(const std::vector<std::string>& cmdline, const NTupleEvent* event)
{
  static bool firstinit = true;
  if (firstinit) {
    firstinit = false;

    FILE* hammer_file = fopen(cmdline[0].c_str(), "r");
    if (not hammer_file) {
      return false;
    }
    Py_Initialize();

    const int argc = cmdline.size();
    char** g_argv = new char*[argc];
    for (int i = 0; i < argc; i++) {
      g_argv[i] = const_cast<char*>(cmdline[i].c_str());
    }
    for (int i = 0; i < argc; i++) {
      printf("%d : '%s'\n", i, g_argv[i]);
    }
    PySys_SetArgv(argc, g_argv);
    PyRun_SimpleString("import ROOT\nimport hammer\n");
  }

  PyObject* py_argv = PyList_New(cmdline.size());
  for (int i = 0; i < cmdline.size(); i++) {
    PyList_SetItem(py_argv, i, PyString_FromString(cmdline[i].c_str()));
  }
  PyObject* py_funcargs = PyTuple_Pack(1, py_argv);
  PyObject* py_func = PyObject_GetAttrString(PyImport_AddModule("hammer"), "main");
  PyObject* py_selector = PyObject_CallObject(py_func, py_funcargs);

  PyRun_SimpleString("hammer.selector_list[-1].SlaveBegin()\n");

  PyObject* py_Global = PyModule_GetDict(PyImport_AddModule("__main__"));
  PyObject* py_Local = PyModule_GetDict(PyImport_AddModule("hammer"));
  PyObject* py_value = PyRun_String("selector_list[-1].this_to_int()", Py_eval_input, py_Global, py_Local);
  SelectorCommon* current = reinterpret_cast<SelectorCommon*>(PyInt_AsLong(py_value));
  selector_list.push_back(current);
  Py_DECREF(py_value);

  current->Init(event);

  return true;
}

bool RootAnalysis::Analyse(const NTupleEvent& event)
{
  for (unsigned n = 0; n < selector_list.size(); n++) {
    SelectorCommon* selector = selector_list[n];

    selector->event_trials = event.trials;
    selector->Process();
  }
  return true;
}

bool RootAnalysis::Finish()
{
  PyRun_SimpleString(
"for selector in hammer.selector_list:\n"
"  selector.SlaveTerminate()\n"
"  selector.stat_report()\n");
//   Py_Finalize();  // will break Sherpa NTuple output if uncommented
  if (g_argv) {
    delete[] g_argv;
  }
  return true;
}
