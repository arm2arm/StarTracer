// CAT2MYSQL.cpp : Defines the entry point for the console application.
//

#include <functional>
#include <algorithm>
#include <string>
#include <utility>
#include <map>

#include <iostream>
#include <sstream> 
#include <string>  
#include <ctime>
#include <ctype.h>
#include <list>
// Needed for the ifstream class
#include <fstream>
// Needed for the setw stream manipulator
#include <iomanip>
using namespace std;
const string mySQLInsert("INSERT INTO `startracer_box64lg`.`st_halos` ( \
`stepnumber` ,\
`id` ,\
`grpid` ,\
`flag` ,\
`xCenter` ,\
`yCenter` ,\
`zCenter` ,\
`radius` ,\
`Npart` ,\
`mass` ,\
`particleCount` ,\
`rmax` ,\
`Rx` ,\
`Ry` ,\
`Rz` ,\
`Ngas` ,\
`Ndm` ,\
`Nst`) VALUES ( ");
void usage(const char*execname)
	{
	cout<<"\n****************************************************"<<endl;
	cout<<"A.Khalatyan"<<endl;
	cout<<"cat2mysql v0.1 2008, Marseille\n"<<endl;
	cout<<"Convertor from AFOF catalog to mysql \"INSERT\" format file."<<endl;
	cout<<"This output can be directly used by mysql<outfile.mysql command to fill tables.\n"<<endl;
	cout<<"Usage:"<<endl;
	cout<<execname<<" [FILEIN.grp.DAT] \n"<<endl;
	cout<<"program output will be ASCII file in mysql format: [FILEIN.grp.DAT].mysql"<<endl;
	cout<<"****************************************************\n"<<endl;

	}
bool StringToInt(const string &s, int &i)
{
  istringstream myStream(s);
  
  if (myStream>>i)
    return true;
  else
    return false;
}

bool GetSnap( string snap, int &isnap)
	{
		string  dig_snap;
	basic_string <char>::size_type indexCh;
    static const basic_string <char>::size_type npos = -1;

	indexCh=snap.find_last_of('/');
	if((indexCh)!=npos)
		snap.assign(snap,int(indexCh+1), int(snap.size()-1));
	
	while(snap.find_first_of("_")!=-1)
		{
			indexCh=snap.find_first_of("_");
			dig_snap=dig_snap.assign(snap,int(indexCh+1), int(3));
			if (StringToInt(dig_snap, isnap))
				return true;
			else
				snap.assign(snap,int(indexCh+1), int(snap.size()-1));
 
		}
	return false;
	}

bool Catalog2MySQL(string file)
	{
	int Nhalo=0;
	float x[3], R,Rxyz[3], flag=1;
	int ID;
	int grpID,Ntotal, npart[3];
	int isnap=0;
	ifstream file_to_read(file.data());
	const int max_num_of_char_in_a_line = 512,// Maximum number of characters expected in a single line in the header
		num_of_header_lines = 1; // Number of header files to skip
	cout<<"Reading catalog file: "<<file;
	if(file_to_read.fail())
		return false;
	///////////////////////////////
	string snap(file);
	if(!GetSnap(snap, isnap))return false;
	///////////////////////////////
	file_to_read.ignore(1, '#');
	file_to_read>>Nhalo;
	file_to_read.ignore(max_num_of_char_in_a_line, '\n');
	int tint=0;
	cout<<endl;
	snap=file+string(".mysql");
	FILE *fout=fopen(snap.c_str(), "w");
	for( int ih=0;ih<Nhalo;ih++)
		{

		file_to_read>>ID;
		file_to_read>>grpID;
		file_to_read>>x[0]>>x[1]>>x[2];

		file_to_read>>R;
		file_to_read>>Ntotal;
		file_to_read>>npart[0]>>npart[1]>>npart[2];
		file_to_read>>Rxyz[0]>>Rxyz[1]>>Rxyz[2];
		file_to_read.ignore(max_num_of_char_in_a_line, '\n');
		if(!file_to_read.good()){file_to_read.close();return false;}
		
		fprintf(fout,"%s %d,%d,%d,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f, %d,%d,%d);\n",mySQLInsert.c_str(), 
			isnap,ID, grpID,flag, x[0], x[1], x[2], 
			R, Ntotal,float(Ntotal),Ntotal,R,
			Rxyz[0], Rxyz[1], Rxyz[2], 
			npart[0], npart[1], npart[2]);
		}
	fclose(fout);
	///////////////////////////////

	file_to_read.close();
	return true;
	}
int main(int argc, char* argv[])
	{
	if(argc!=2)
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

		}
	else
		{
		string fname=string(argv[1]);
		if(!Catalog2MySQL(fname))
			{
			cout<<"Error fail to read and convvert...exiting..."<<endl;
			exit(EXIT_FAILURE);
			}
		cout<<"done"<<endl;
		}

	return EXIT_SUCCESS;
	}

