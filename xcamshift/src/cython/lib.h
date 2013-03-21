// -*- c++ -*-
#ifndef FWK_LIB_H
#define FWK_LIB_H 1

struct _object;
typedef _object PyObject;

class IAlg 
{
 public:
  IAlg();
  virtual ~IAlg();
  void doRun();
  virtual void run() =0;
  void doJump(int i);
  virtual void jump(int i);
};

IAlg *create_alg(int algid);

class CyAlgBase : public IAlg
{
public:
  PyObject *m_obj;

  CyAlgBase(PyObject *obj);
  virtual ~CyAlgBase();
  virtual void run();
  virtual void jump(int i);
};

#endif // !FWK_LIB_H
