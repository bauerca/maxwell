#include "MxUtil.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "hdf5.h"

#include "MxShape.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h"

#include "MxMap.hpp"
#include "MxCrsMatrix.hpp"

typedef std::pair<std::string, std::string> MxStrPair;

char listSpacers[] = {',', ' ', '[', ']', '{', '}'};


int MxUtil::pow(int base, int exp) {
  int res = 1;
  for (int i = 0; i < exp; ++i)
    res *= base;
  return res;
}

template<typename T, int DIM>
MxDimVector<double, DIM> MxUtil::ndRootFind(const std::vector<const T *> & funcObjs, MXDIMVEC guess, double tol, int maxiter) {
  int jacRows = funcObjs.size();
  if (jacRows != DIM)
    std::cout << "MxUtil::ndRootFind(...): over or under-determined system.\n";

  MxDimMatrix<double, DIM> jac;
  MxDimVector<double, DIM> f, x, dx;

  x = guess;
  for (int i = 0; i < maxiter; ++i) {

    // find func vals and build jacobian at trial point
    for (int j = 0; j < jacRows; ++j) {
      f[j] = funcObjs[j]->func(x);
      std::cout << "f[" << j << "] = " << f[j] << "\n";
      jac.setRow(j, funcObjs[j]->gradFunc(x));
    }
    std::cout << "Jacobian:\n"; jac.print();
    std::cout << "Inverse Jacobian:\n"; jac.inv().print();

    if (f.oneNorm() < tol) {
      std::cout << "|f| = " << f.oneNorm() << "\n";
      return x;
    }

    jac.solve(dx, f);

    if (dx.oneNorm() < tol) {
      std::cout << "|dx| = " << dx.oneNorm() << "\n";
      return x;
    }

    x += dx;
  }

  std::cout << "MxUtil::ndRootFind(...): maxiter exceeded. Returning root with\n"
            << "  |dx| = " << dx.oneNorm() << "\n"
            << "   |f| = " << f.oneNorm() << "\n";
  return x;
}

template MxVecD3 MxUtil::ndRootFind<MxShape<3>, 3>(const std::vector<const MxShape<3> *> & funcObjs, MxVecD3 guess, double tol, int maxiter);
template MxVecD2 MxUtil::ndRootFind<MxShape<2>, 2>(const std::vector<const MxShape<2> *> & funcObjs, MxVecD2 guess, double tol, int maxiter);

template<>
double MxUtil::Strings::strToType<double>(const std::string & str) {
  return atof(str.c_str());
}

template<>
std::string MxUtil::Strings::strToType<std::string>(const std::string & str) {
  return str;
}

template<>
int MxUtil::Strings::strToType<int>(const std::string & str) {
  return atoi(str.c_str());
}

// default is for T a number 
template<typename T>
std::vector<T> MxUtil::Strings::strToStdVector(const std::string & str) {
  //std::cout << "strToStdVector str = " << str << "\n";
  std::vector<T> res;
  std::string::const_iterator iter;

  bool prevWasSpacer = true;
  bool curIsSpacer = true;

  std::vector<char> vecListSpacers(listSpacers, listSpacers + 4);

  std::string substr;
  for (iter = str.begin(); iter != str.end(); ++iter) {
    curIsSpacer = MxUtil::inStdVector(*iter, vecListSpacers);
    if (!curIsSpacer) {
      substr.append(1, *iter);
      if (iter == str.end() - 1) {
        //std::cout << "  " << substr << "\n";
        res.push_back(MxUtil::Strings::strToType<T>(substr));
      }
    }
    else if (!prevWasSpacer) {
      //std::cout << "  " << substr << "\n";
      res.push_back(MxUtil::Strings::strToType<T>(substr));
      substr.clear();
    }
    prevWasSpacer = curIsSpacer;
  }

  //MxUtil::printStdVector(res);

  return res;
}

template std::vector<double> MxUtil::Strings::strToStdVector<double>(const std::string & str);
template std::vector<int> MxUtil::Strings::strToStdVector<int>(const std::string & str);
template std::vector<std::string> MxUtil::Strings::strToStdVector<std::string>(const std::string & str);

#if 0
template<>
std::vector<std::string> MxUtil::Strings::strToStdVector<std::string>(const std::string & str) {
  std::vector<std::string> res;
  std::string::const_iterator iter;

  bool prevWasSpacer = true;
  bool curIsSpacer = true;

  std::vector<char> vecListSpacers(listSpacers, listSpacers + 4);

  std::string substr;
  for (iter = str.begin(); iter != str.end(); ++iter) {
    curIsSpacer = MxUtil::inStdVector(*iter, vecListSpacers);
    if (!curIsSpacer) {
      substr.append(1, *iter);
      if (iter == str.end() - 1)
        res.push_back(substr); // implicit conversion from double to T
    }
    else if (!prevWasSpacer) {
      res.push_back(substr); // implicit conversion from double to T
      substr.clear();
    }
    prevWasSpacer = curIsSpacer;
  }

  return res;
}
#endif

std::string MxUtil::Strings::stripLeadingSpaces(std::string const & str) {
  size_t pos = str.find_first_not_of(" ");
  return str.substr(pos);
}

std::string MxUtil::Strings::stripTrailingSpaces(std::string const & str) {
  size_t pos = str.find_last_not_of(" ");
  return str.substr(0, pos + 1);
}

std::string MxUtil::Strings::stripLTSpaces(std::string const & str) {
  std::string tmp = MxUtil::Strings::stripLeadingSpaces(str);
  return MxUtil::Strings::stripTrailingSpaces(tmp);
}

std::string MxUtil::Strings::lhsEquation(std::string const & str) {
  return MxUtil::Strings::splitEquation(str).first;
}

std::pair<std::string, std::string> MxUtil::Strings::splitEquation(std::string const & str) {
  std::pair<std::string, std::string> res;
  std::string left, right;

  size_t pos = str.find('=');

  left = str.substr(0, pos);
  right = str.substr(pos + 1);
  //std::cout << "left: " << left << "    right: " << right << "\n";

  // now get rid of leading/trailing spaces in each expression
  res.first = MxUtil::Strings::stripLTSpaces(left);
  res.second = MxUtil::Strings::stripLTSpaces(right);

  return res;
}

std::vector<std::string> MxUtil::Strings::lines(std::string const & str) {
  std::vector<std::string> res;
  std::string line;
  std::string::const_iterator iter;
  for (iter = str.begin(); iter != str.end(); ++iter) {
    if (*iter == '\n') {
      res.push_back(line);
      line.clear();
    }
    else
      line.append(1, *iter);
  }
  return res;
}

std::vector<std::string> MxUtil::Strings::split(std::string const & str, std::string sep) {
  std::vector<std::string> res;
  size_t pos1 = 0, pos2 = 0;
  while (pos2 != std::string::npos) {
    pos2 = str.find(sep, pos1);
    res.push_back(str.substr(pos1, pos2 - pos1));
    pos1 = pos2 + sep.size();
  }
  return res;
}

// this function gathers all lines within an xml object block that are not contained within
// a nested xml object block. Ignores blank lines.
std::vector<std::string> MxUtil::XML::nodeLines(const ::Teuchos::XMLObject & node) {
  std::vector<std::string> res;
  std::vector<std::string> allLines = MxUtil::Strings::lines(node.toString());
  size_t numTotLines = allLines.size();
  std::string line;
  int numNodesOpen = 0;
  bool blank;

  for (int i = 0; i < numTotLines; ++i) {
    line = allLines[i];
    blank = line.find_first_not_of(" ") == std::string::npos;

    if (line.find("</") != std::string::npos)
      numNodesOpen--;
    else if (line.find("<") != std::string::npos)
      numNodesOpen++;
    else if (numNodesOpen == 1 and not blank)
      res.push_back(line);
    else
      continue;
  }

  //MxUtil::printStdVector(res);

  return res;
}


std::string MxUtil::XML::getAttr(std::string name, const ::Teuchos::XMLObject & node, std::string dfault) {
  std::string resB, resT;

  if (node.isEmpty())
    return dfault;

  resB = MxUtil::XML::getBodyAttr(name, node);
  resT = MxUtil::XML::getTagAttr(name, node);

  if (resB != "") return resB;
  else if (resT != "") return resT;
  else return dfault;
}

// If attr is not found, this function returns an error and exits.
std::string MxUtil::XML::getAttr(std::string name, const ::Teuchos::XMLObject & node) {
  std::string s = MxUtil::XML::getAttr(name, node, "");
  if (s == "") {
    //return an error
    std::cout << "MxUtil::XML::getAttr(...): attribute, '" << name << "', is required to exist in the following XML block:\n";
    node.print(std::cout, 1);
    exit(0);
  } 
  else
    return s;
}

std::string MxUtil::XML::getBodyAttr(std::string name, const ::Teuchos::XMLObject & node) {
  std::vector<std::string> lines = MxUtil::XML::nodeLines(node);
  MxStrPair eq;

  std::vector<std::string>::const_iterator iter;
  for (iter = lines.begin(); iter != lines.end(); ++iter) {
    eq = MxUtil::Strings::splitEquation(*iter);
    //std::cout << "eq.first: " << eq.first << "   eq.second: " << eq.second << "\n";
    if (eq.first == name)
      return eq.second;
  }

  // couldn't find 'name' so we return an empty string
  return "";
}

// this one's easy with Trilinos
std::string MxUtil::XML::getTagAttr(std::string name, const ::Teuchos::XMLObject & node) {
  return node.getWithDefault(name, std::string(""));
}

template<typename Scalar>
void MxUtil::Trilinos::massiveCrsMultiply(
std::vector<RCP<MxCrsMatrix<Scalar> > > const & mats,
std::vector<bool> const & transposes,
//std::vector<char> const & transposes,
const RCP<MxCrsMatrix<Scalar> > & res, bool fc) {

  size_t numMats = mats.size();
  int pid = res->getRangeMap()->getComm()->myPID();

  RCP<MxCrsMatrix<Scalar> > m1 = mats[0];
  //char trans1 = transposes[0];
  bool trans1 = transposes[0];

  RCP<const MxCrsMatrix<Scalar> > m2; // *m2 is always one of the given matrices
  //char trans2;
  bool trans2;

  RCP<MxCrsMatrix<Scalar> > tmpres; // structure of tmpres is determined from m2/m1

  if (pid == 0) std::cout << "Multiplying " << numMats << " matrices:\n";
  for (size_t i = 1; i < numMats; i++) {
    m2 = mats[i];
    trans2 = transposes[i];
    if (i == numMats - 1)
      tmpres = res;
    else
      //tmpres = rcp(new MxCrsMatrix<Scalar>(
      //  trans2 ? m2->getDomainMap() : m2->getRangeMap(),
      //  trans1 ? m1->getRangeMap() : m1->getDomainMap()));
      tmpres = rcp(new MxCrsMatrix<Scalar>(
        trans2 ? m2->getDomainMap() : m2->getRangeMap()));
      //tmpres = rcp(new MxCrsMatrix<Scalar>(trans2 == 'T' ?
      //  m2->getDomainMap() : m2->getRangeMap()));

    if (pid == 0) std::cout << "  multiplying...\n";
    if (i == numMats - 1)
      MxCrsMatrix<Scalar>::matrixMatrixMultiply(
          *m2, trans2, *m1, trans1, *tmpres, fc);
      //MxCrsMatrix<Scalar>::matrixMatrixMultiply(
      //    *m2, trans2 == 'T' ? true : false,
      //    *m1, trans1 == 'T' ? true : false, *tmpres, fc);
    else
      MxCrsMatrix<Scalar>::matrixMatrixMultiply(
          *m2, trans2, *m1, trans1, *tmpres, true);
      //MxCrsMatrix<Scalar>::matrixMatrixMultiply(
      //    *m2, trans2 == 'T' ? true : false,
      //    *m1, trans1 == 'T' ? true : false, *tmpres, true);

    m1 = tmpres;

    //trans1 = 'F';
    trans1 = false;
  }
}

template void MxUtil::Trilinos::massiveCrsMultiply<double>(
std::vector<RCP<MxCrsMatrix<double> > > const &,
std::vector<bool> const &,
const RCP<MxCrsMatrix<double> > &, bool);
template void MxUtil::Trilinos::massiveCrsMultiply<MxComplex>(
std::vector<RCP<MxCrsMatrix<MxComplex> > > const &,
std::vector<bool> const &,
const RCP<MxCrsMatrix<MxComplex> > &, bool);

#if 0
template<typename Scalar>
void MxUtil::Trilinos::massiveMultiply(
std::vector<RCP<MxCrsMatrix<Scalar> > > mats,
std::vector<char> transposes,
RCP<MxCrsMatrix<Scalar> > & res, bool fc) {

  size_t numMats = mats.size();
  int pid = res.Comm().MyPID();

  Epetra_CrsMatrix * m1 = new Epetra_CrsMatrix(*mats[0]);
  char trans1 = transposes[0];
  const Epetra_CrsMatrix * m2; // *m2 is always one of the given matrices
  char trans2;
  Epetra_CrsMatrix * tmpres; // structure of tmpres is determined from m2/m1

  if (pid == 0) std::cout << "Multiplying " << numMats << " matrices:\n";
  for (size_t i = 1; i < numMats; i++) {
    m2 = mats[i];
    trans2 = transposes[i];
    if (i == numMats - 1)
      tmpres = &res;
    else
      tmpres = new Epetra_CrsMatrix(Copy, trans2 == 'T' ? m2->DomainMap() : m2->RangeMap(), 0);

    //if (trans1 && trans2) {
    //  if (!m2->RangeMap().SameAs(m1->DomainMap())) {
    //    std::cout << "Matrix " << i-1
    if (pid == 0) std::cout << "  multiplying...\n";
    if (i == numMats - 1)
      EpetraExt::MatrixMatrix::Multiply(*m2, trans2 == 'T' ? true : false, *m1, trans1 == 'T' ? true : false, *tmpres, fc);
    else
      EpetraExt::MatrixMatrix::Multiply(*m2, trans2 == 'T' ? true : false, *m1, trans1 == 'T' ? true : false, *tmpres, true);

    delete m1; //because m1 will soon point to result of above line
    m1 = tmpres; //here

    trans1 = 'F';
  }
}
#endif

template<typename Scalar>
RCP<MxCrsMatrix<Scalar> > MxUtil::Trilinos::identity(RCP<MxMap> map) {
  RCP<MxCrsMatrix<Scalar> > res;
  res = rcp(new MxCrsMatrix<Scalar>(map, map)); // static profile 
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  size_t n = map->getNodeNumIndices();
  MxIndex gid;

  for (size_t i = 0; i < n; i++) {
    gid = map->getGlobalIndex(i);
    res->insertRowValues(gid, 1, &gid, &one);
  }

  res->fillComplete(map, map);
  return res;
}

template RCP<MxCrsMatrix<double> > MxUtil::Trilinos::identity<double>(
RCP<MxMap> map);

template RCP<MxCrsMatrix<MxComplex> > MxUtil::Trilinos::identity<MxComplex>(
RCP<MxMap> map);



void MxUtil::Epetra::stripZeros(Epetra_CrsMatrix & m) {

  if (!m.Filled()) {
    std::cout << "MxUtil::Epetra::stripZeros: input matrix has to be FillCompleted.";
    throw 1;
  }

  Epetra_Map dMap(m.DomainMap()), rMap(m.RangeMap());

  Epetra_CrsMatrix copy(Copy, m.RowMap(), m.ColMap(), 0);

  //std::cout << "Global inds m?: " << m.IndicesAreGlobal() << "\n";
  //std::cout << "Global inds copy?: " << copy.IndicesAreGlobal() << "\n";

  //strip zeros from resulting matrix (why doesn't Trilinos do this?) to conserve space
  int numRowInds = m.RowMap().NumMyElements();
  int numEntries;
  int * indices;
  double * values;
  for (int i = 0; i < numRowInds; i++) {
    m.ExtractMyRowView(i, numEntries, values, indices);
    for (int j = 0; j < numEntries; j++) {
      if (fabs(values[j]) > 1.e-12)
        copy.InsertMyValues(i, 1, &values[j], &indices[j]);
    }
  }
  copy.FillComplete(dMap, rMap);

  m = copy;
}



void MxUtil::Epetra::removeConstField(Epetra_MultiVector & x) {
  //x.Comm().Barrier();
  int nVecs = x.NumVectors();

  Epetra_Vector ones(x.Map());
  ones.PutScalar(1);

  //double dotProds[nVecs];
  //double norms[nVecs];
  double scale, dotProd;

  //std::cout << "norm = " << norm << ", dotProd = " << dotProd << "\n";
  for (int i = 0; i < nVecs; ++i) {
    ones.Dot(ones, &scale);
    x(i)->Dot(ones, &dotProd);
    //std::cout << "const vec part of input vector: " << dotProd << "\n";
    x(i)->Update(-dotProd / scale, ones, 1.);
    x(i)->Dot(ones, &dotProd);
    //std::cout << "const vec part after projection: " << dotProd << "\n";
  }
}

void MxUtil::Epetra::printNorms(std::string name, Epetra_MultiVector const & x) {
  std::vector<double> norm2(x.NumVectors());
  x.Norm2(&norm2[0]);

  std::cout << name << " norms: ";
  for (int i = 0; i < x.NumVectors(); i++)
    std::cout << norm2[i] << ", ";
  std::cout << "\n";
}

void MxUtil::HDF5::saveArray(double * data, int length, const char * name) {
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  hsize_t dimsf = length;

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  // for hdf5 >= 1.8
#if H5_VERS_MINOR == 8
  hid_t llist_id = H5Pcreate(H5P_LINK_CREATE);
  hid_t alist_id = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dset_id = H5Dcreate2(file_id, "data", H5T_NATIVE_DOUBLE, filespace_id, llist_id, H5P_DEFAULT, alist_id);
  H5Pclose(llist_id);
  H5Pclose(alist_id);
  // for hdf5 1.6
#else
  hid_t dset_id = H5Dcreate(file_id, "data", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT);
#endif

  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace_id);
  H5Fclose(file_id);
}



