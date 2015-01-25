// StarTracer.cpp : Defines the entry point for the console application.
//
#pragma warning(disable:4996)

#include <functional>
#include <algorithm>
#include <cstring>
#include <string>
#include <utility>
#include <map>

#include <iostream>
#include <sstream> 
#include <stdlib.h>
#ifdef WIN32
#include <conio.h>
#endif
#include <stdio.h>
#include <ctime>
#include <ctype.h>
#include <list>
#include <cmath>
// Needed for the ifstream class
#include <fstream>
// Needed for the setw stream manipulator
#include <iomanip>

#include "tree.hh"
#include "tree_util.hh"

using namespace std;
#define DEBUG_ME true
#define RadFrac  0.6

//////////////////////
bool StringToInt(const string &s, int &i)
{
  istringstream myStream(s);
  
  if (myStream>>i)
    return true;
  else
    return false;
}
bool GetSnapPath( string snap, string &path)
	{
	basic_string <char>::size_type indexCh;
	static const basic_string <char>::size_type npos = 0;

	indexCh=snap.find_last_of('/');
	if((indexCh)!=npos)
		path.assign(snap,0, int(indexCh));
	else
		path.assign("./");
	return true;
	}
bool GetSnap( string &snap, int &isnap)
	{
	string  dig_snap;
	basic_string <char>::size_type indexCh;
	static const basic_string <char>::size_type npos = 0;

	indexCh=snap.find_last_of('/');
	if((indexCh)!=npos)
		snap.assign(snap,int(indexCh+1), int(snap.size()-1));
	//cout<<"Snap: "<<snap<<endl;
	while(snap.find_last_of("_")!=0)
		{
		indexCh=snap.find_last_of("_");
		dig_snap=dig_snap.assign(snap,int(indexCh+1), int(4));
		//	cout<<dig_snap<<endl;
		if (StringToInt(dig_snap, isnap))
			{
			snap=dig_snap;
			return true;
			}
		else
			snap.assign(snap,int(indexCh+1), int(snap.size()-1));

		}
	return false;
	}
//////////////////////

class CTimer
	{
	public:
		void start(){time (&m_start);};
		void stop(bool print_flag=false){
			time (&m_end);
			m_dif = difftime (m_end,m_start);
			if(print_flag)
				print();
			}
		void print(void ){cout<<m_dif<<" secs"<<endl;};
		time_t m_start;
		time_t m_end;
		double m_dif;
	};
CTimer timer;
// set_intersection example
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;
using std::ifstream;

#define  int4bytes int
/* PRINT_ELEMENTS()
* - prints optional C-string optcstr followed by
* - all elements of the collection coll
* - separated by spaces
*/
template <class T>
inline void PRINT_ELEMENTS (const T& coll, const int nel=10,const char* optcstr="")
	{
	typename T::const_iterator pos;
	int i=0;
	std::cout << optcstr;
	int np=min(int(coll.size()), int(nel));
	for (pos=coll.begin(); i<np; ++pos, ++i) {
		std::cout << *pos << endl;
		}
	std::cout << std::endl;
	}

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
	{
	size_t nread;

	if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
		{
		printf("I/O error (fread) %d !\n",(unsigned int )nread);
		exit(3);
		}
	return nread;
	}

struct strParticleData
	{
	 float Pos[3];
	};

class CReaders
	{
	public:
bool good(){return m_isgood;};	
	protected:
		bool m_isgood;
};

class CData :public CReaders
	{
	public:
		CData():swp_flag(false),find_flag(false)
			{};
		virtual ~CData() {};
		virtual bool ReadData(string file){return true;};
		int GetBlk(ifstream *file_to_read, int *blk)
			{

			file_to_read->rdbuf()->sgetn(reinterpret_cast<char *>(blk), (std::streamsize)sizeof(int));
			bool bad=file_to_read->bad();
			return *blk;
			};
		void GetBlkName(ifstream *file_to_read, char *name)
			{

			file_to_read->rdbuf()->sgetn(name,sizeof(char)*4);
			GetBlk(file_to_read, &blk);
			};

		void swap_Nbyte(char *data,int n,int m)
			{
			int i,j;
			char old_data[16];

			if(swp_flag)
				{
				for(j=0;j<n;j++)
					{
					memcpy(&old_data[0],&data[j*m],m);
					for(i=0;i<m;i++)
						{
						data[j*m+i]=old_data[m-i-1];
						}
					}
				}
			};
		size_t my_fread(void *ptr, size_t size, size_t nmemb, ifstream *file_to_read)
			{
			size_t nread;

			if((nread = file_to_read->rdbuf()->sgetn((char*)ptr,size*nmemb)) != size*nmemb)
				{
				if(!find_flag)
					{
					cerr<<"I/O error (fread) "<<nread<<endl;
					exit(3);
					}else
					{
					cerr<<"Cannot find block: ";
					return 0;
						}
				}
			return nread;
			};

		bool swp_flag;
		bool find_flag;
		int blk;
	};

class CGadget : public  CData
	{
	public:
		struct io_header
			{
			int npart[6];
			double mass[6];
			double time;
			double redshift;
			int flag_sfr;
			int flag_feedback;
			int npartTotal[6];
			int flag_cooling;
			int num_files;
			double BoxSize;
			double Omega0;
			double OmegaLambda;
			double HubbleParam;
			int flag_multiphase;
			int flag_stellarage;
			int flag_sfrhistogram;
			int flag_metals;
			int flag_decouple;
			int flag_effmodel;
			char fill[72];		/* fills to 256 Bytes */
			}myhead;

		CGadget(string file):P(NULL), ID(NULL)
			{
			m_isgood=ReadData(file);
			};
		~CGadget(){
			if(ID!=NULL)
				delete ID;
			if(P!=NULL)
				delete P;
			};
		bool ReadData(string file);
		void GetHeader(ifstream *file_to_read);
		bool good(){return m_isgood;};	
		int find_block(ifstream *fd,char *label)
			{
			int4bytes blocksize=0, blksize=0, ret=0;
			char blocklabel[5]={"    "};
			find_flag=true;
			fd->seekg(ios::beg);
			//  printf("Finding: %s\n",label);
			while(!fd->eof() && (blocksize == 0))//&& strcmp(blocklabel, "Z   ")!=0 )
				{
				//      cout<<"SKIP1"<<endl;
				GetBlk(fd, &blksize);
				//      cout<<"SKIP2"<<endl;
				//				if(blksize == 134217728)
					{				
					swap_Nbyte((char*)&blksize,1,4);
					}
					if(blksize != 8)
						{
						if(ret>0)
							{
							printf("incorrect format (blksize=%d)!\n",blksize);
							exit(1);
							}else break;
						}
					else
						{

						ret=my_fread(blocklabel, 4*sizeof(char), 1, fd);
						if(ret>0)
							ret=my_fread(&blocksize, sizeof(int4bytes), 1, fd);
						swap_Nbyte((char*)&blocksize,1,4);
						if( DEBUG_ME)
							{
							if(ret>0)
								printf("Found Block <%s> with %d bytes\n",blocklabel,blocksize);
							else
								printf(" <%s> \n",label);

							}

						GetBlk(fd, &blksize);
						if(strcmp(label,blocklabel)!=0)
							{ 	
							fd->seekg(blocksize,ios_base::cur);
							blocksize=0;
							}
						}
				}
			find_flag=false;
			return(blocksize-8);
			}
			void SeekToType(ifstream *file_to_read,int type, int onesize)
			{
			 long nskip=0;
			 for(int ip=0;ip<type;ip++)
				nskip+=this->myhead.npart[ip];
			 file_to_read->seekg((nskip)*onesize, ios_base::cur);
			}

		int *ID;
		strParticleData *P;
		bool m_isgood;

	};
void CGadget::GetHeader(ifstream *fd){
	GetBlk(fd, &blk);
	my_fread((void*)myhead.npart,6*sizeof(int), 1, fd);             swap_Nbyte((char*)myhead.npart,6,4);
	my_fread((void*)myhead.mass,6*sizeof(double), 1, fd);           swap_Nbyte((char*)myhead.mass,6,8);
	my_fread((void*)&myhead.time,sizeof(double), 1, fd);            swap_Nbyte((char*)&myhead.time,1,8);
	my_fread((void*)&myhead.redshift,sizeof(double), 1, fd);        swap_Nbyte((char*)&myhead.redshift,1,8);
	my_fread((void*)&myhead.flag_sfr,sizeof(int), 1, fd);           swap_Nbyte((char*)&myhead.flag_sfr,1,4);
	my_fread((void*)&myhead.flag_feedback,sizeof(int), 1, fd);      swap_Nbyte((char*)&myhead.flag_feedback,1,4);
	my_fread((void*)myhead.npartTotal,6*sizeof(int), 1, fd);        swap_Nbyte((char*)myhead.npartTotal,6,4);
	my_fread((void*)&myhead.flag_cooling,sizeof(int), 1, fd);       swap_Nbyte((char*)&myhead.flag_cooling,1,4);
	my_fread((void*)&myhead.num_files,sizeof(int), 1, fd);          swap_Nbyte((char*)&myhead.num_files,1,4);
	my_fread((void*)&myhead.BoxSize,sizeof(double), 1, fd);         swap_Nbyte((char*)&myhead.BoxSize,1,8);
	my_fread((void*)&myhead.Omega0,sizeof(double), 1, fd);          swap_Nbyte((char*)&myhead.Omega0,1,8);
	my_fread((void*)&myhead.OmegaLambda,sizeof(double), 1, fd);     swap_Nbyte((char*)&myhead.OmegaLambda,1,8);
	my_fread((void*)&myhead.HubbleParam,sizeof(double), 1, fd);     swap_Nbyte((char*)&myhead.HubbleParam,1,8);
	my_fread((void*)&myhead.flag_multiphase,sizeof(int), 1, fd);    swap_Nbyte((char*)&myhead.flag_multiphase,1,4);
	my_fread((void*)&myhead.flag_stellarage,sizeof(int), 1, fd);    swap_Nbyte((char*)&myhead.flag_stellarage,1,4);
	my_fread((void*)&myhead.flag_sfrhistogram,sizeof(int), 1, fd);  swap_Nbyte((char*)&myhead.flag_sfrhistogram,1,4);
	my_fread((void*)&myhead.flag_metals,sizeof(int), 1, fd);        swap_Nbyte((char*)&myhead.flag_metals,1,4);
	my_fread((void*)&myhead.flag_decouple,sizeof(int), 1, fd);      swap_Nbyte((char*)&myhead.flag_decouple,1,4);
	my_fread((void*)&myhead.flag_effmodel,sizeof(int), 1, fd);      swap_Nbyte((char*)&myhead.flag_effmodel,1,4);
	my_fread((void*)myhead.fill,72*sizeof(char), 1, fd);
	GetBlk(fd, &blk);
cout << "=====================" << endl;
    for (unsigned it = 0; it < 6; it++) {

        printf("N[%d]=%0.9d\tMass[%d]=%g\n", it, myhead.npart[it], it,
                myhead.mass[it]);
    };
    cout << "=====================" << endl;
	};


bool CGadget::ReadData(string file)
	{
	ifstream file_to_read(file.data(),  ios::in|ios::binary);
	bool isgood=true, is_multi=false;
	char name[5];
	unsigned long Ntotal;
	vector<string> vmmfile;
	memset(name, 0,sizeof(name));
	if(file_to_read.bad())
	       isgood=false;
	else
	  {
	    file_to_read.close();
	    string mmfile=file+".0";
	    file_to_read.open(mmfile.data(), ios::in|ios::binary);
	    
	    if(!file_to_read.bad()){
	      isgood=true;
	      is_multi=true;
	      cout<<"Found multiple file format"<<endl;

	    }else 
	      return false;
	  }

	GetBlk(&file_to_read, &blk);
	if(blk != 8)
		{
		swp_flag=true;
		swap_Nbyte((char*)&blk,1,4);
		if(blk!=8)return false;
		}

	GetBlkName(&file_to_read, name);

	GetBlk(&file_to_read, &blk);
	GetHeader(&file_to_read);
	if(is_multi){
	  Ntotal=myhead.npartTotal[4];
	  char buffer[4];
	  for(int i=0;i<myhead.num_files;i++){
	    sprintf(buffer, ".%d", i);
	    vmmfile.push_back(file + string(buffer));
	  }
	}
	else
	  Ntotal=myhead.npart[4];
	
	ID=new int [Ntotal];
	P=new strParticleData [Ntotal];
	size_t icurr_pos=0;
	for(int fi=0;fi<myhead.num_files;fi++){

	  cout<<"Reading: "<<vmmfile[fi]<<endl;
	  int nid=0, id=0, sizeall, size, nskip=0;
	  file_to_read.close();
	  file_to_read.open(vmmfile[fi].data(), ios::in|ios::binary);
	  
	  find_block(&file_to_read, (char *)"HEAD");
	  GetHeader(&file_to_read);
	  cout<<this->myhead.npart[4]<<endl;

	  GetBlk(&file_to_read, &nid);
	  sizeall=find_block(&file_to_read, (char *)"ID  ");
	  size=this->myhead.npart[4]*sizeof(int);
	  
	  if(sizeall !=size){
	    SeekToType(&file_to_read,4, sizeof(int));
	  }
	  
	  GetBlk(&file_to_read, &blk);
	  my_fread(ID+icurr_pos, size, 1, &file_to_read);
	  swap_Nbyte((char*)(ID+icurr_pos),this->myhead.npart[4],4);
	  
	  
	  sizeall=find_block(&file_to_read, (char *)"POS ");
	  GetBlk(&file_to_read, &blk);
	  size=this->myhead.npart[4]*sizeof(float);
	  
	  
	  SeekToType(&file_to_read,4, sizeof(strParticleData));
	  my_fread(P+icurr_pos, size, 3, &file_to_read);
	  swap_Nbyte((char*)(P+icurr_pos),this->myhead.npart[4]*3,4);
	  icurr_pos+=this->myhead.npart[4];
	  
	}

return true;
}

class CGalaxy;

typedef pair <vector<CGalaxy*>::iterator, float> Map_Int_Flt_Pair;

//////////////////////////////////////////
class CGalaxy{
public:
	CGalaxy(int idgrp=0):grpID(idgrp){
		memset(x,0,sizeof(x));
		//		memset(Rxyz,0,sizeof(Rxyz));
		memset(npart,0,sizeof(npart));
		Ntotal=0;R=0;ID=0;m_Rxx=0;
		};

	~CGalaxy(){};
	void SetIds(int nids){id.resize(nids,0);};
	void SortIds()
		{
		if(id_rxx.size()>0)
			sort(id_rxx.begin(), id_rxx.end());
			if(id.size()>0)
				{
					sort(id.begin(), id.end());
					int nold=id.size();
					id.erase (unique(id.begin(),
                        id.end()),
	                 id.end());
					if(nold!=id.size())
						{
						cout << "Corrupted ID list for IDGrp="<<this->ID<<" "<<this->grpID<<" found "<<
	nold-id.size()<<" duplicates out of "<<nold<<" in ID list: "<< endl;
						//exit(0);
						}
					
					
				}
		};
	/////////////////////////////
	inline bool operator == (const CGalaxy &b) const
		{
		return (b.grpID==grpID);
		}
	/////////////////////////////
	//vector<CIntersect> child_list;
	vector<Map_Int_Flt_Pair> child_list;
	vector<int> id;
	vector<int> id_rxx;
	float x[3];
	float R;
	float m_Rxx;
	int ID;
	int grpID;
	int npart[3];
	//int Rxyz[3];
	int Ntotal;
protected:
	};

//////////////////////////////////////
//////////////////////////////////////////
// Comparators for CGalaxy objects
// Ascending date sorting function
class cmpADMassSort
	{
	public:
		cmpADMassSort(bool mode=true):m_mode(mode){};
		bool operator()( CGalaxy*const& rpStart,  CGalaxy*const& rpEnd)
			{
			if(m_mode)
				return (rpStart->Ntotal) < (rpEnd->Ntotal);
			else
				return (rpStart->Ntotal) > (rpEnd->Ntotal);
			}
	private: 
		bool m_mode;
	};

// Descending date sorting function
struct cmpDescendingMassSort
	{
	bool operator()(CGalaxy*& rpStart, CGalaxy*& rpEnd)
		{
		return rpStart->Ntotal > rpEnd->Ntotal;
		}
	};
//////////////////////////////////////

struct HalogrpId: public std::binary_function< CGalaxy, int, bool > {
	bool operator () ( const CGalaxy &gal, const int &grpID ) const {
		return gal.grpID == grpID;
		}
	};


// this is our function object class.
// it holds a reference to a vector of grpIDS.
class CgrpIDScomp : public std::binary_function<unsigned int,unsigned int,bool> 
	{
	// the vector reference
	const vector<unsigned int>	&m_grpIDS;

	public:
		// constructor which takes a reference to a vector of CgrpIDS
		CgrpIDScomp( const vector<unsigned int> & CgrpIDS ) : m_grpIDS(CgrpIDS) {}

		// comparison operator. this will be called by std::sort with two
		// integers. we'll perform a comparison operation on the inputs and
		// return the result of the comparison.
		bool operator()(int a, int b) const
			{
			// a typical comparison operator compares the two input parameters.
			// this one uses the two input parameters as indexes into the m_grpIDS vector.
//			return (m_grpIDS.at(a).ID) < (m_grpIDS.at(b).ID);

			return (m_grpIDS.at(a)) < (m_grpIDS.at(b));
			}
	};

class CFOFCatalog : public CReaders
	{
	public:
		CFOFCatalog(){};
		CFOFCatalog(string file);
		~CFOFCatalog(){};
		bool ReadCatalog(string file);
		bool ReadCatalogIDs(string file);
		bool ReadCatalogIDsBin(string file);
		bool InitIDs( CGadget *gad)
			{
			Redshift=(float)gad->myhead.redshift;
			m_gadobj=gad;
			return InitIDs( gad->ID);
			}
		bool InitIDs( int *ID)
			{
			cout<<"Init IDS ...";
			bool ret=false;
			string file=m_afofcat, filebin;
			file.resize(m_afofcat.size()-4);
			IDs=ID;
			filebin=file;
			filebin+=string(".bin");
			cout<<".RB.";
			ret=ReadCatalogIDsBin(filebin);
			if(!ret)
				{
				cout<<".RA.";
				ret=ReadCatalogIDs(file);
				}
			cout<<"ok"<<endl;
			return ret;
			}
		void FillIDs();

		vector<CGalaxy*> halos;
		vector<unsigned int> grpIDs;
		int *IDs;
		CGadget* m_gadobj;
		float Redshift;

		vector<int> indexVector;// used for sorting 
		vector<int>::iterator it;// used for accessing to sorted elements
		string m_afofcat;
		
	};

CFOFCatalog::CFOFCatalog(string file):m_afofcat(file)
	{
	cerr<<"Reading catalog file :"<<file<<endl;
	if(!(m_isgood=ReadCatalog(file)))
		{
		cerr<<"Cannot open catalog file :"<<file<<endl;
		cerr<<"Exiting..."<<endl;
		}


	}
bool CFOFCatalog::ReadCatalogIDsBin(string filename)
	{
	unsigned int nid=0, i=0;
	
	ifstream ifile(filename.c_str(),ios::in|ios::binary);
	if(	ifile.fail())
		{
		string mes=string("Can not open file: " );
		mes+= filename;
		cout<<mes.c_str()<<endl;
		return false;
		}

	ifile.read((char*)(&nid), sizeof(unsigned int));
	grpIDs.resize(nid);
	cout<<"..";
	ifile.read((char*)(&grpIDs[0]), sizeof(unsigned int)*nid);
	ifile.close();
	indexVector.clear();
	for(i=0;i<nid;i++)
		{
		// initialize indexVector
		if(grpIDs[i] >0)
			{
			indexVector.push_back(i);
			}
		}

	timer.start();
	cout<<".S.";
	sort(indexVector.begin(), 
		indexVector.end(), 
		CgrpIDScomp(grpIDs)); 
	timer.stop(DEBUG_ME);
	FillIDs();
	return true;
	};


class CPosDist2comp : public std::binary_function<int,int,bool> 
	{
	const vector<float>	&m_grpIDS;
	public:
		CPosDist2comp( const vector<float> & CgrpIDS ) : m_grpIDS(CgrpIDS) {}
		bool operator()(int a, int b) const
			{
			return (m_grpIDS.at(a)) < (m_grpIDS.at(b));
			}
	};

void CFOFCatalog::FillIDs()
	{
	unsigned int iindex=0, ind=0, ii=0;
	static bool first_flag=true;
	vector<CGalaxy*>::iterator it;// used for accessing to sorted elements
	cout<<".FID. "<<indexVector.size()<<"  ";
#ifdef DumpIDs
	FILE *testf;
	
	if(first_flag)
		{
		testf=fopen("C:/mingw/1.0/home/arm2arm/DATA/pos.txt", "w");
		first_flag=false;
		}else
	  testf=fopen("C:/mingw/1.0/home/arm2arm/DATA/pos.txt", "a");
	
#endif
		vector<float> dist2;
		vector<int> indexOfDist2;
		float d2;
	for(iindex=0;iindex<indexVector.size();iindex++)
		{
		int value=grpIDs[indexVector[iindex]];

		for(it=halos.begin();it<halos.end();it++)
			{
			if((*it)->grpID == value)
				{
				(*it)->id.resize((*it)->Ntotal);
				dist2.clear();
				indexOfDist2.clear();
				for(ind=0;ind<(unsigned int)((*it)->Ntotal);ind++)
					{
					(*it)->id[ind]=IDs[indexVector[iindex+ind]];
					strParticleData *p=&(this->m_gadobj->P[indexVector[iindex+ind]]);
					indexOfDist2.push_back(ind);
					d2=0;
					for(int i=0;i<3;i++)				
						d2+= ((*p).Pos[i]-(*it)->x[i])*((*p).Pos[i]-(*it)->x[i]);//Sum(x-xc)^2
					dist2.push_back(d2);

					}
				sort(indexOfDist2.begin(), indexOfDist2.end(), CPosDist2comp(dist2));
				int iRpos=int(dist2.size()*RadFrac);
				(*it)->m_Rxx=sqrt(dist2.at(iRpos));
#ifdef DumpIDs
				if((*it)->ID==1)
					{
					fprintf(testf,"#%d  \n", iRpos);
				for(ind=0;ind<(unsigned int)((*it)->Ntotal);ind++)
					{
					    strParticleData *p=&(this->m_gadobj->P[ indexVector[iindex+indexOfDist2[ind]] ] );
						fprintf(testf,"%f %f %f \n", (*p).Pos[0],(*p).Pos[1],(*p).Pos[2]);
					
					};
					}
#endif
				// Copy vector truncated by Rxx 
				(*it)->id_rxx.resize(iRpos);
				for(ind=0;ind<(unsigned int)(iRpos);ind++)
					{
						(*it)->id_rxx[ind]=(*it)->id[indexOfDist2[ind]];
					}
				//sort and check  
				(*it)->SortIds();
				iindex+=(*it)->Ntotal-1;
				//cout<<ii++<<"("<<(*it)->Ntotal<<")"<<"  ";
				}

			}
		}
#ifdef DumpIDs	
	fprintf(testf,"\n\n");
	fclose(testf);
#endif
	
	}
bool CFOFCatalog::ReadCatalogIDs(string file)
	{

	ifstream file_to_read(file.data());
	//	const int max_num_of_char_in_a_line = 512,// Maximum number of characters expected in a single line in the header
	//	num_of_header_lines = 1; // Number of header files to skip

	timer.start();
	if(DEBUG_ME)cout<<"Reading ID file: "<<file;
	if(file_to_read.bad())
		return false;
	int nid=0, id=0,i;
	
	file_to_read>>nid;
	///////////////////////////////
	grpIDs.resize(nid);
	///////////////////////////////
	for(i=0;i<nid;i++)
		{
		file_to_read>>grpIDs[i];
		// initialize indexVector
		if(grpIDs[i] >0)
			{
			indexVector.push_back(i);
			}
		}
	///////////////////////////////
	file_to_read.close();
	if(DEBUG_ME)
		cout<<"\ntotal "<<indexVector.size()<<"  and greather than zero: "<<float(indexVector.size())/float(nid)*100<<"perc... ok"<<endl;
	///////////////////////////////
	timer.stop(DEBUG_ME);
	// Do index sort by grpIDs
	if(DEBUG_ME)
		{
		int iindex=0;
		for(iindex=0;iindex<10;iindex++)
			cout<<indexVector[iindex]<<" "<<grpIDs[indexVector[iindex]]<<endl;


		cout<<"make heap"<<endl;
		timer.start();	
		//make_heap(indexVector.begin(), indexVector.end(),CgrpIDScomp(grpIDs));
		timer.stop(DEBUG_ME);		
		cout<<"sorting heap"<<endl;
		}
	timer.start();

	sort(indexVector.begin(), 
		indexVector.end(), 
		CgrpIDScomp(grpIDs)); 
	timer.stop(DEBUG_ME);
	FillIDs();
	return true;


	}

bool CFOFCatalog::ReadCatalog(string file)
	{
	ifstream file_to_read(file.data());
	const int max_num_of_char_in_a_line = 512,// Maximum number of characters expected in a single line in the header
		num_of_header_lines = 1; // Number of header files to skip
	if(DEBUG_ME)
		cout<<"Reading catalog file: "<<file;
	if(file_to_read.bad())
		return false;

	///////////////////////////////
	int Nhalo=0;
	file_to_read.ignore(1, '#');
	file_to_read>>Nhalo;
	file_to_read.ignore(max_num_of_char_in_a_line, '\n');
	CGalaxy* g;
	int tint=0;
	for( int ih=0;ih<Nhalo;ih++)
		{
		g=new CGalaxy;

		file_to_read>>g->ID;
		file_to_read>>g->grpID;
		file_to_read>>g->x[0]>>g->x[1]>>g->x[2];

		file_to_read>>g->R;
		file_to_read>>g->Ntotal;
		file_to_read>>g->npart[0]>>g->npart[1]>>g->npart[2];
		//		   file_to_read>>g->Rxyz[0]>>g->Rxyz[1]>>g->Rxyz[2];
		file_to_read.ignore(max_num_of_char_in_a_line, '\n');
		if(!file_to_read.good())return false;
		halos.push_back(g);
		}
	///////////////////////////////

	file_to_read.close();
	if(DEBUG_ME)
		cout<<"... ok"<<endl;

	if(DEBUG_ME)
		{
		cout<<"before sort"<<endl;
		for(int i=0;i<10;i++)
			cout<<(halos[i])->Ntotal<<endl;
		}
	if(DEBUG_ME)
		cout<<"Sorting catalog...";

	sort(halos.begin(), halos.end(),cmpADMassSort(false));

	if(DEBUG_ME)
		cout<<"...ok"<<endl;

	if(DEBUG_ME)
		{
		cout<<"after sort"<<endl;	
		for(int i=0;i<10;i++)
			cout<<(halos[i])->Ntotal<<endl;
		}

	return true;
	}

class CGenTree
	{
	public:
		CGenTree(){delete m_cat1; delete m_cat2;};
		CGenTree(CFOFCatalog *cat1, CFOFCatalog *cat2):m_cat1(cat1),m_cat2(cat2) {};
		void MakeTree(void );
		void MakeTreeExt(void);
		void DumpTree(void);
		void DumpProgenitorsTree(void);
		void OpenFile();
	protected:
		FILE *ftree;
		CFOFCatalog *m_cat1;
		CFOFCatalog *m_cat2;
		int isnap1, isnap2;
		tree<CGalaxy*> m_tree;
	};

const char *MYSQL_UPDATE_processed_steps=
"INSERT INTO `startracer_box64lg`.`processed_steps` (\
`startStep` ,\
`finalStep` ,\
flaggedDone) \
VALUES (\
%d, %d, 1);\n";
const char *MYSQL_UPDATE_tree_processed_steps=
"INSERT INTO `startracer_box64lg`.`processed_tree_steps` (\
`stepNumber`) \
VALUES (\
%d);\n";

const char *MYSQL_UPDATE_available_steps=
"INSERT INTO `startracer_box64lg`.`available_steps` (\
`stepNumber` ,\
flagged, Redshift) \
VALUES (\
%d, 1, %f);\n";

const char *MYSQL_Insert_progenitors=
"INSERT INTO `startracer_box64lg`.`halos_progenitors` (\
`startStep` ,\
`finalStep` ,\
`startHalo` ,\
`finalHalo` ,\
`progenitor_particle_count` ,\
`progenitor_mass` ,\
`inclusion_rate` ,\
`reverse_inclusion_rate`\
)\
VALUES (\
%d, %d, %d, %d, %d, %f, %f, %f\
);\n";

const char *MYSQL_Insert_progenitors_tree=
"INSERT INTO `startracer_box64lg`.`progenitors_tree` (\
`parentStepNumber` ,\
`parentHaloId` ,\
`childStepNumber` ,\
`childHaloId` ,\
`orderCol` ,\
`childMass` )\
VALUES (\
%d, %d, %d, %d, %d, %f);\n";

void CGenTree::OpenFile()
	{
	static bool first_time=true;
	string snap1(m_cat1->m_afofcat);
	string snap2(m_cat2->m_afofcat);
	string path, fnametree(m_cat1->m_afofcat);
	if(!GetSnap( snap1, isnap1)){cout<<"cannot get snap number form the filename..."<<endl;exit(13);};
	if(!GetSnap( snap2, isnap2)){cout<<"cannot get snap number form the filename..."<<endl;exit(13);};
	GetSnapPath( fnametree,  path);
	fnametree=path+string("/tree_")+snap1+string("_")+string(snap2)+string(".mysql");	
	if(first_time)
		{
		ftree=fopen(fnametree.c_str(), "w");
		first_time=false;
		}else
		ftree=fopen(fnametree.c_str(), "a");
	cout<<endl;
	cout<<fnametree<<"  .Wr."<<flush;
	}
void CGenTree::DumpTree(void)
	{
	OpenFile();
	
	fprintf(ftree,MYSQL_UPDATE_processed_steps, isnap1, isnap2);
	fprintf(ftree,MYSQL_UPDATE_tree_processed_steps, isnap1);
	fprintf(ftree,MYSQL_UPDATE_available_steps, isnap1, m_cat1->Redshift);

	unsigned int total=m_cat1->halos.size();
	for(unsigned int i=0;i<total;i++)
		{
		if((i%50==0))cout<<"."<<flush;
		unsigned int total_child=m_cat1->halos[i]->child_list.size();
		//if(!total_child)cout<<"ERROR..... "<<i<<" "<<m_cat1->halos[i]->ID<<" "<<m_cat1->halos[i]->Ntotal<<endl;;
		for(unsigned int ich=0;ich<total_child;ich++)
			{
			vector<CGalaxy*>::iterator it_ch=m_cat1->halos[i]->child_list[ich].first, 
				it=(m_cat1->halos.begin()+i);
			float inclusion_rate=m_cat1->halos[i]->child_list[ich].second;
			fprintf(ftree,MYSQL_Insert_progenitors, isnap1, isnap2, (*it)->ID,
				(*it_ch)->ID, (*it_ch)->Ntotal,(float)(*it_ch)->Ntotal, inclusion_rate,inclusion_rate);
				//m_cat1->halos[i]->child_list[ich].second<<") "<<flush;
			}
		
		}
	cout<<".ok"<<endl;
	fclose(ftree);
	};		
void CGenTree::DumpProgenitorsTree(void)
	{
	OpenFile();
	
	unsigned int total=m_cat1->halos.size();
	for(unsigned int i=0;i<total;i++)
		{
		if((i%50==0))cout<<"."<<flush;
		unsigned int total_child=m_cat1->halos[i]->child_list.size();
		//if(!total_child)cout<<"ERROR..... "<<i<<" "<<m_cat1->halos[i]->ID<<" "<<m_cat1->halos[i]->Ntotal<<endl;;
		for(unsigned int ich=0;ich<total_child;ich++)
			{
			vector<CGalaxy*>::iterator it_ch=m_cat1->halos[i]->child_list[ich].first, 
				it=(m_cat1->halos.begin()+i);
			float inclusion_rate=m_cat1->halos[i]->child_list[ich].second;
			fprintf(ftree,MYSQL_Insert_progenitors_tree, isnap1,(*it)->ID, isnap2, 
				(*it_ch)->ID, (*it_ch)->Ntotal,(float)(*it_ch)->Ntotal);
			}
		
		}
	cout<<".ok"<<endl;
	fclose(ftree);
	};		

void CGenTree::MakeTree(void)
	{

	vector<int> res;
	list<int> L;
	vector<int> arr;
	vector<int>::iterator arr_end;
	list<vector<CGalaxy*>::iterator> gal2_list;
	unsigned int iCount=0;
	timer.start();
	cout<<endl;
	cout<<".MkTree.";
	for(unsigned int i=0;i<m_cat1->halos.size();i++)
		{
		if((i%50)==0)
			cout<<"."<<flush;

		if(m_cat1->halos[i]->Ntotal > 0)
			{

			for(unsigned int j=0;j<m_cat2->halos.size();j++)
				{
				iCount=0;
				arr.clear();
				arr.resize(min(m_cat1->halos[i]->Ntotal,m_cat2->halos[j]->Ntotal));
				arr_end=set_intersection(m_cat1->halos[i]->id.begin(),
					m_cat1->halos[i]->id.end(),
					m_cat2->halos[j]->id.begin(),
					m_cat2->halos[j]->id.end(),
					arr.begin());
				
				iCount=(unsigned int)(arr_end-arr.begin());
				if(iCount>0)
					{
					//cout<<"  N1="<<m_cat1->halos[i]->Ntotal<<"  N2="<<m_cat2->halos[j]->Ntotal;
					//cout<<"intersected FW= "<< iCount<<endl;
					m_cat1->halos[i]->child_list.push_back(
						Map_Int_Flt_Pair((m_cat2->halos.begin()+j),
						iCount/float(m_cat1->halos[i]->Ntotal)*100
						)
						);//end of insert
					}


				}
			}
		
		}
	
	cout<< ".done   time: ";
	timer.stop(true);
	
	}
void CGenTree::MakeTreeExt(void)
	{

	vector<int> res;
	list<int> L;
	vector<int> arr;
	vector<int>::iterator arr_end;
	list<vector<CGalaxy*>::iterator> gal2_list;
	unsigned int iCount=0;
	timer.start();
	cout<<endl;
	cout<<".MkTree.";
	for(unsigned int i=0;i<m_cat1->halos.size();i++)
		{
		if((i%50)==0)
			cout<<"."<<flush;

		if(m_cat1->halos[i]->Ntotal > 0)
			{

			for(unsigned int j=0;j<m_cat2->halos.size();j++)
				{
				iCount=0;
				arr.clear();
				arr.resize(min(m_cat1->halos[i]->Ntotal,m_cat2->halos[j]->Ntotal));
				arr_end=set_intersection(m_cat1->halos[i]->id.begin(),
					m_cat1->halos[i]->id.end(),
					m_cat2->halos[j]->id.begin(),
					m_cat2->halos[j]->id.end(),
					arr.begin());
				
				iCount=(unsigned int)(arr_end-arr.begin());
				if(iCount>0)
					{
					//cout<<"  N1="<<m_cat1->halos[i]->Ntotal<<"  N2="<<m_cat2->halos[j]->Ntotal;
					//cout<<"intersected FW= "<< iCount<<endl;
					m_cat1->halos[i]->child_list.push_back(
						Map_Int_Flt_Pair((m_cat2->halos.begin()+j),
						iCount/float(m_cat1->halos[i]->Ntotal)*100
						)
						);//end of insert
					}


				}
			}
		
		}
	
	cout<< ".done   time: ";
	timer.stop(true);
	
	}
void usage(const char*execname)
	{
	cout<<"\n*****************************************************"<<endl;
	cout<<"A.Khalatyan"<<endl;
	cout<<"StarTracer v0.1 2008, Marseille\n"<<endl;
	cout<<"Creates mysql query file to fill\nAFOF groups host-progenotor relations."<<endl;
	cout<<"This tables can be used directly by \nthe Tree-web interface.\n"<<endl;
	cout<<"Usage:"<<endl;
	cout<<execname<<" [snappath_snap1_]\\ \n [snappath_snap1_]\\ \n [grppath_snap1_.grp.DAT]\\ \n [grppath_snap2_.grp.DAT] \n"<<endl;
	cout<<"program output will be ASCII file: \n\t  [grppath]/tree_snap1_snap2.mysql"<<endl;
	cout<<"Comments: to fill table use: \n\t  mysql -uMySqlUsername < \\ \n\t\t[grppath]/tree_snap1_snap2.mysql"<<endl;
	cout<<"*******************************************************\n"<<endl;
	}

int main(int argc, char **argv)
	{
	if(argc!=5)
		{
		//Under linux directory separator given by "/"
		char* pch=strrchr(argv[0], '/');
		//Under windows directory separator given by "\"
		if(pch==NULL)
			pch=strrchr(argv[0], '\\');
		if(pch!=NULL)
			{	string execname=string(pch+1);
		usage(execname.c_str());
			}else
				usage(argv[0]);

		}else
		{

		string SnapName1=string(argv[1]);
		string SnapName2=string(argv[2]);
		string CatName1=string(argv[3]);
		string CatName2=string(argv[4]);
		
		/*		int isnap;
		GetSnap( SnapName1, isnap);
		cout<<isnap<<endl;exit(0);
		*/
		
		CGadget gad(SnapName1);
		CGadget gad2(SnapName2);
		CFOFCatalog cat1(CatName1);
		CFOFCatalog cat2(CatName2);
		bool goodness=gad.good() && gad2.good()&& cat1.good() && cat2.good();
		if(goodness)
			{
			cat1.InitIDs(&gad);
			cat2.InitIDs(&gad2);

			CGenTree gentree(&cat1, &cat2);
			gentree.MakeTree();
			gentree.DumpTree();
			gentree.DumpProgenitorsTree();
			cout<<endl<<"all done."<<endl;
#ifdef WIN32
		int ch = _getch();
#endif

			return EXIT_SUCCESS;
			}
		}

		cout<<"\nsomething is odd...please check input files."<<endl;
		return EXIT_FAILURE;

	}
