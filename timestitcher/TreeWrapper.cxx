#include <iostream>

#include "TreeWrapper.h"

using namespace std;

namespace R3BCalifaTimestitcher
{

   TreeWrapper::TreeWrapper(TTree *tree, uint32_t id) : tree(tree), idx(0), id(id)
   {
      this->nEntries = tree->GetEntries();
   }

   uint64_t TreeWrapper::getTS()
   {
      return *ptrTS;
   }

   bool TreeWrapper::next()
   {
//      cerr << "TreeWrapper::next()" << endl;
//      Long64_t treeIdx;
      while(idx < nEntries)
      {
//         cerr << " idx = " << idx << endl;
//         treeIdx = tree->LoadTree(idx++);   
//         cerr << " treeIdx = " << treeIdx << endl;
         tree->GetEntry(idx++);
         if(isGood())
            return true;
      }
      return false;
   }

   uint32_t TreeWrapper::getId()
   {
      return id;
   }

}
