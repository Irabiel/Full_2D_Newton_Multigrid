#include "tools.h"

#include <sys/stat.h>
#include <string>

int mkDir(const char *path)
{
#ifdef _WIN32
    return ::_mkdir(path);
#else
//#if _POSIX_C_SOURCE
//    return ::mkdir(path);
//#else
    return ::mkdir(path, 0755);
//#endif
#endif
};

int GetFileLen(FILE* myFile)
{
    fseek (myFile, 0, SEEK_END);
    int size = ftell(myFile);
    fseek(myFile, 0, SEEK_SET);
    return size;
};

char* GetNextString(char*& buffer)
{
    char* out = buffer;
    if (!*buffer) return NULL; // return on empty string
    while(! (*buffer == 0x0A || *buffer == 0x0D || *buffer == 0x00) ) // 0x0A and 0x0D
        buffer++; // skip forward until we find the start of the next line (10/13/0)
    if (*buffer) *buffer++ = 0; // if we ended on 10/13 end the string and move to the next char
    if(*buffer == 0x0A) buffer++;  // on windows skip the 10 after the 13
    
    return out;
};

param ReadParameters(char* fname)
{
    FILE* FID = fopen(fname, "r");
    if (FID == NULL) {
        fprintf(stderr, "Can't open parameter file.\n");
        exit(1);
    }
    
    char* data_string;
    char var_name[100];
    char var_value[100];
    
    int fileLen = GetFileLen(FID);
    char* buffer = (char*) malloc(fileLen+1);
    fread(buffer, fileLen, 1, FID);
    buffer[fileLen] = 0;
    
    param p;
    
    while((data_string = GetNextString(buffer)))
    {
        sscanf(data_string, "%s %s", var_name, var_value);
        
        if (strcmp(var_name,"BoxX")==0)
            p.BoxX = atof(var_value);
        else if (strcmp(var_name,"BoxY")==0)
            p.BoxY = atof(var_value);
        else if (strcmp(var_name,"Levels")==0)
            p.levels = atof(var_value);
		else if (strcmp(var_name,"Cycles")==0)
            p.Cycles = atof(var_value);
		else if (strcmp(var_name,"A0")==0)
            p.A0 = atof(var_value);
		else if (strcmp(var_name,"B0")==0)
            p.B0 = atof(var_value);
        else if (strcmp(var_name,"DirName")==0)
        {
            strcpy(p.DirName, var_value);
            mkDir(p.DirName);
        }
        else
        {
            printf("Unknown parameter: %s \n", var_name);
            /*fflush(stdout);
             assert(false);
             exit(-1);*/
        }
    }
    fclose(FID);
    
    return p;
};
