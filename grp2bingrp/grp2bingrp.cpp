// grp2bingrp.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <sstream> 
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <ctime>
using namespace std;
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
void usage(const char*execname)
	{
	cout<<"\n****************************************************"<<endl;
	cout<<"A.Khalatyan"<<endl;
	cout<<"grp2bingrp v0.1 2008, Marseille\n"<<endl;
	cout<<"Convertor for ASCII AFOF grp file to binary format.\nThis output can be directly used by StarTracer.\n"<<endl;
	cout<<"Usage:"<<endl;
	cout<<execname<<" [FILEIN.grp] \n"<<endl;
	cout<<"program output will be binary file: [FILEIN.grp].bin"<<endl;
	cout<<"****************************************************\n"<<endl;

	}
class CDataConvertor{
public:
	CDataConvertor(string infn):infilename(infn){
		outfilename=infilename;
		outfilename+=string(".bin");
		data=NULL;
		memblock=NULL;
		};
	~CDataConvertor(){
		if (memblock!=NULL)
			delete [] memblock;
		if (data!=NULL)
			delete [] data;
		}
	void DoConvert()
		{
		cout<<".R."<<flush();
		ReadBuffer(infilename);
		cout<<".P."<<flush();
		ParseBuffer();
		cout<<".W."<<flush();
		WriteBuffer(outfilename);
		cout<<"..."<<flush();
		};
	string infilename;
	string outfilename;
	char * memblock;
	unsigned int *data; 
	unsigned int nelem;
	long size;
	void error(const char *msg)
		{
		cout<<"\nError:\n"<<msg<<endl;
		cout<<"exiting ....\n"<<endl;
		exit(EXIT_FAILURE);
		};
	long GetFileSize(ifstream *infile)
		{
		long sizeoffile, begin, end;
		infile->seekg (0, ios::beg);
		begin = infile->tellg();
		infile->seekg (0, ios::end);
		end = infile->tellg();
		sizeoffile=end-begin;
		string msg=string("The file ");
		msg+=string(" is EMPTY!!!");
		if( sizeoffile<1)error(msg.c_str());
		infile->seekg (0, ios::beg);
		return sizeoffile;
		};
	void ParseBuffer(bool timing_flag=false)
		{
		long n=0;
		stringstream stream(memblock);
		stream>>nelem;
		if(timing_flag)
			{
			timer.start();
			cout<<nelem<<"  ";
			}
		data=new unsigned int[nelem];
		if(data==NULL)
			error("Can not allocate memory for data in function: ParseBuffer");

		for(unsigned int i=0;i<nelem;i++)
			{
			stream>>data[i];
			}
		if(timing_flag)
			{
			timer.stop(true);
			}
		};
	void ReadBuffer(string infilename)
		{
		ifstream infile(infilename.c_str(),ios::in|ios::binary);
		if(	infile.fail())
			{
			string mes=string("Can not find file: " );
			mes+= infilename;

			error(mes.c_str());
			}

		size=GetFileSize(&infile);
		memblock = new char [size];
		infile.seekg (0, ios::beg);
		infile.read (memblock, size);
		infile.close();
		};
	void WriteBuffer(string filename)
		{
		long size=0;
		ofstream ofile(filename.c_str(),ios::out|ios::binary);
		if(	ofile.fail())
			{
			string mes=string("Can not open file: " );
			mes+= filename;
			error(mes.c_str());
			}

		ofile.write((char*)(&nelem), sizeof(unsigned int));
		ofile.write((char*)(&data[0]), sizeof(unsigned int)*nelem);
		ofile.close();
		};

	};
int main(int argc, char* argv[])
	{
	if(argc!=2)
		{
		//Under linux directory separator given by "/"
		char* pch=strrchr(argv[0], '/');
		//Under windows directory separator given by "\"
		if(pch==NULL)
			pch=strrchr(argv[0], '\\');
		string execname=string(pch+1);
		usage(execname.c_str());
		}
	else
		{
		string fname=string(argv[1]);
		CDataConvertor *convertor=new CDataConvertor(fname);
		cout<<"converting ascii file: ";
		cout<<convertor->infilename;
		cout<<" to ";
		cout<<convertor->outfilename;
		cout<<"...";
		convertor->DoConvert();
		delete convertor;
		cout<<"done"<<endl;
		}

	return EXIT_SUCCESS;
	}

