
#include <Python.h>
#include <cstdio>
#include <iostream>

#include "SelectorCommon.h"

#include "root-analysis.h"

using namespace RootAnalysis;

static char** g_argv = 0;

static std::vector<SelectorCommon*> selector_list;

int RootAnalysis::Init(const std::vector<std::string>& cmdline)
{
  static bool firstinit = true;
  if (firstinit) {
    firstinit = false;

    FILE* hammer_file = fopen(cmdline[0].c_str(), "r");
    if (not hammer_file) {
      throw;
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

  PyRun_SimpleString("hammer.selector_list[-1].SlaveBegin(0)\n");

  PyObject* py_Global = PyModule_GetDict(PyImport_AddModule("__main__"));
  PyObject* py_Local = PyModule_GetDict(PyImport_AddModule("hammer"));
  PyObject* py_value = PyRun_String("selector_list[-1].this_to_int()", Py_eval_input, py_Global, py_Local);
  SelectorCommon* current = reinterpret_cast<SelectorCommon*>(PyInt_AsLong(py_value));
  selector_list.push_back(current);
  Py_DECREF(py_value);

  return 0;
}

int RootAnalysis::Analyse(const NTupleEvent& event)
{
  for (unsigned n = 0; n < selector_list.size(); n++) {
    SelectorCommon* selector = selector_list[n];

    selector->event_trials = event.trials;

    selector->ntuple_id = event.id;
    selector->ntuple_nparticle = event.nparticle;
    assert(event.nparticle <= MAXNPARTICLE);
    for (int i = 0; i < event.nparticle; i++) {
      selector->ntuple_kf[i] = event.kf[i];
      selector->ntuple_px[i] = event.px[i];
      selector->ntuple_py[i] = event.py[i];
      selector->ntuple_pz[i] = event.pz[i];
      selector->ntuple_E[i] = event.E[i];
    }
    selector->ntuple_alphas = event.alphas;
    selector->ntuple_weight = event.weight;
    selector->ntuple_me_wgt = event.me_wgt;
    selector->ntuple_x1 = event.x1;
    selector->ntuple_x2 = event.x2;
    selector->ntuple_x1p = event.x1p;
    selector->ntuple_x2p = event.x2p;
    selector->ntuple_id1 = event.id1;
    selector->ntuple_id2 = event.id2;
    selector->ntuple_fac_scale = event.fac_scale;
    selector->ntuple_ren_scale = event.ren_scale;
    selector->ntuple_nuwgt = event.nuwgt;
    for (int i = 0; i < event.nuwgt; i++) {
      selector->ntuple_usr_wgts[i] = event.usr_wgts[i];
    }
    selector->ntuple_alphaspower = event.alphaspower;
    selector->ntuple_part[0] = event.part[0];
    selector->Process(-1);
  }
}

int RootAnalysis::Finish()
{
  PyRun_SimpleString(
"for selector in hammer.selector_list:\n"
"  selector.SlaveTerminate()\n"
"  selector.stat_report()\n");
//   Py_Finalize();  // will break Sherpa NTuple output if uncommented
  if (g_argv) {
    delete[] g_argv;
  }
  return 0;
}
