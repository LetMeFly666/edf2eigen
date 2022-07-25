/*
***************************************************************************
*
* Author: Teunis van Beelen
*
* Copyright (C) 2007 - 2019 Teunis van Beelen
*
* Email: teuniz@protonmail.com
*
***************************************************************************
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
***************************************************************************
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>

#if defined(__APPLE__) || defined(__MACH__) || defined(__APPLE_CC__)

#define fopeno fopen

#else

#define fseeko fseek
#define ftello ftell
#define fopeno fopen

#endif


#define ELECTRODE_TAG "[ELECTRODE]"
#define ELECTRODE_UNTAG "["
#define ELECTRODE_NAME_MAXLEN 256

#define ANNOT_TRACKSIZE 54


int total_elapsed_time;

char labels[256][17];


int check_device(char *);

int convert_nk2edf(FILE *, FILE *, FILE *,  int, int, int, char *, int);

void latin1_to_utf8(char *, int);

void latin1_to_ascii(char *, int);

int read_21e_file(char *);



int main(int argc, char *argv[])
{
    FILE *inputfile=NULL,
          *outputfile=NULL,
           *logfile=NULL,
            *pntfile=NULL;

    const char *fileName;

    int i, j,
        pathlen,
        fname_len,
        error,
        ctl_block_cnt,
        datablock_cnt,
        total_blocks,
        edfplus,
        n_logs=0,
        n_sublogs=0,
        total_logs=0,
        n_logblocks=0,
        ctlblock_address,
        wfmblock_address,
        logblock_address,
        read_subevents=0;

    char path[512],
         logfilepath[512],
         pntfilepath[512],
         *log_buf=NULL,
         *sublog_buf=NULL;

    total_elapsed_time = 0;

    setlocale(LC_NUMERIC, "C");

    if((argc!=2)&&(argc!=3)) {
        printf("\nNihon Kohden to EDF(+) converter. ver. 1.5\n"
               "Copyright 2007 - 2019 Teunis van Beelen\n"
               "Email: teuniz@protonmail.com\n"
               "This software is licensed under the GNU GENERAL PUBLIC LICENSE Version 3.\n\n"
               "Usage: nk2edf [-no-annotations] <filename>\n\n"
               "normal use: nk2edf <filename>\n"
               "Three files are needed with the extions .eeg, .pnt and .log\n"
               "this will create an EDF+ file including annotations.\n\n"
               "A *.21E file will be used (if present) to read in alternative electrode names.\n\n"
               "If only the .eeg file is available, use: nk2edf -no-annotations <filename>\n"
               "this will create an EDF file without annotations.\n\n");
        return(1);
    }

    if(argc==2) {
        strcpy(path, argv[1]);
        edfplus = 1;
    } else {
        strcpy(path, argv[2]);
        if(strcmp(argv[1], "-no-annotations")) {
            printf("\nNihon Kohden to EDF(+) converter. ver. 1.4\n"
                   "Copyright 2007 - 2019 Teunis van Beelen\n"
                   "Email: teuniz@protonmail.com\n"
                   "This software is licensed under the GNU GENERAL PUBLIC LICENSE Version 3.\n\n"
                   "Usage: nk2edf [-no-annotations] <filename>\n\n"
                   "normal use: nk2edf <filename>\n"
                   "Three files are needed with the extions .eeg, .pnt and .log\n"
                   "this will create an EDF+ file including annotations.\n\n"
                   "A *.21E file will be used (if present) to read in alternative electrode names.\n\n"
                   "If only the .eeg file is available, use: nk2edf -no-annotations <filename>\n"
                   "this will create an EDF file without annotations.\n\n");
            return(1);
        }
        edfplus = 0;
    }

    pathlen = strlen(path);

    if(pathlen<5) {
        printf("Error, filename must contain at least five characters.\n");
        return(1);
    }

    fname_len = 0;
    for(i=pathlen; i>0; i--) {
        if((path[i-1]=='/')||(path[i-1]=='\\'))  break;
        fname_len++;
    }
    fileName = path + pathlen - fname_len;

    for(i=0; fileName[i]!=0; i++);
    if(i==0) {
        printf("Error, filename must contain at least five characters.\n");
        return(1);
    }

    i -= 4;
    if((strcmp((const char *)fileName + i, ".eeg"))&&(strcmp((const char *)fileName + i, ".EEG"))) {
        printf("Error, filename extension must have the form \".eeg\" or \".EEG\"\n");
        return(1);
    }

    inputfile = fopeno(path, "rb");
    if(inputfile==NULL) {
        printf("Error, can not open file %s for reading\n", path);
        return(1);
    }

    /***************** check if the EEG file is valid ******************************/

    char scratchpad[32];

    rewind(inputfile);
    fread(scratchpad, 16, 1, inputfile);
    scratchpad[16] = 0;
    if(check_device(scratchpad)) {
        printf("Error, deviceblock has unknown signature: \"%s\"\n", scratchpad);
        fclose(inputfile);
        return(1);
    }
    fseeko(inputfile, 0x0081LL, SEEK_SET);
    fread(scratchpad, 16, 1, inputfile);
    scratchpad[16] = 0;
    if(check_device(scratchpad)) {
        printf("Error, controlblock has unknown signature: \"%s\"\n", scratchpad);
        fclose(inputfile);
        return(1);
    }
    fseeko(inputfile, 0x17feLL, SEEK_SET);
    if(fgetc(inputfile)!=0x01) {
        printf("Error, waveformdatablock has wrong signature.\n");
        fclose(inputfile);
        return(1);
    }

    /************************* read logs **********************************************/

    if(edfplus) {
        strncpy(logfilepath, path, 512);
        pathlen = strlen(logfilepath);
        strcpy(logfilepath + pathlen - 3, "log");
        logfile = fopeno(logfilepath, "rb");
        if(logfile==NULL) {
            printf("Can not open file %s for reading,\n"
                   "if there is no .log file you can try to create an EDF file instead of EDF+.\n\n"
                   "Three files are needed with the extions .eeg, .pnt and .log\n"
                   "this will create an EDF+ file including annotations.\n\n"
                   "If only the .eeg file is available, use: nk2edf -no-annotations <filename>\n"
                   "this will create an EDF file.\n\n",
                   logfilepath);
            fclose(inputfile);
            return(1);
        }

        rewind(logfile);
        fread(scratchpad, 16, 1, logfile);
        scratchpad[16] = 0;
        if(check_device(scratchpad)) {
            printf("Error, .log file has unknown signature: \"%s\"\n", scratchpad);
            fclose(logfile);
            fclose(inputfile);
            return(1);
        }

        fseeko(logfile, 0x0091LL, SEEK_SET);
        n_logblocks = fgetc(logfile);
        log_buf = (char *)malloc(n_logblocks * 11521);
        if(log_buf==NULL) {
            printf("Malloc error (logbuf)\n");
            fclose(logfile);
            fclose(inputfile);
            return(1);
        }
        sublog_buf = (char *)malloc(n_logblocks * 11521);
        if(sublog_buf==NULL) {
            printf("Malloc error (sublogbuf)\n");
            fclose(logfile);
            fclose(inputfile);
            free(log_buf);
            return(1);
        }

        read_subevents = 1;

        total_logs = 0;

        for(i=0; i<n_logblocks; i++) {
            fseeko(logfile, 0x0092LL + (i * 20), SEEK_SET);
            if(fread((char *)(&logblock_address), 4, 1, logfile)!=1) {
                printf("Error reading .log file.\n");
                fclose(inputfile);
                fclose(logfile);
                free(log_buf);
                free(sublog_buf);
                return(1);
            }
            fseeko(logfile, logblock_address + 0x0012LL, SEEK_SET);
            n_logs = fgetc(logfile);
            fseeko(logfile, logblock_address + 0x0014LL, SEEK_SET);
            if(fread(log_buf + (total_logs * 45), n_logs * 45, 1, logfile)!=1) {
                printf("Error reading .log file.\n");
                fclose(inputfile);
                fclose(logfile);
                free(log_buf);
                free(sublog_buf);
                return(1);
            }

            if(read_subevents) {
                if(fseeko(logfile, 0x0092LL + ((i + 22) * 20), SEEK_SET)) {
                    read_subevents = 0;
                } else {
                    if(fread((char *)(&logblock_address), 4, 1, logfile)!=1) {
                        read_subevents = 0;
                    } else {
                        if(fseeko(logfile, logblock_address + 0x0012LL, SEEK_SET)) {
                            read_subevents = 0;
                        } else {
                            n_sublogs = fgetc(logfile);
                            if(n_sublogs != n_logs) {
                                read_subevents = 0;
                            } else {
                                if(fseeko(logfile, logblock_address + 0x0014LL, SEEK_SET)) {
                                    read_subevents = 0;
                                } else {
                                    if(fread(sublog_buf + (total_logs * 45), n_sublogs * 45, 1, logfile)!=1) {
                                        read_subevents = 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            total_logs += n_logs;
        }

        for(i=0; i<total_logs; i++) {
            for(j=0; j<20; j++) {
                if(((unsigned char *)log_buf)[(i * 45) + j]<32)  log_buf[(i * 45) + j] = ' ';
            }

            latin1_to_utf8(log_buf + (i * 45), 20);

            if(read_subevents) {
                strncpy(log_buf + (i * 45) + 26, sublog_buf + (i * 45) + 24, 6);
            }
        }

        /************************* check pntfile **********************************************/

        strncpy(pntfilepath, path, 512);
        pathlen = strlen(pntfilepath);
        strcpy(pntfilepath + pathlen - 3, "pnt");
        pntfile = fopeno(pntfilepath, "rb");
        if(pntfile==NULL) {
            printf("Can not open file %s for reading,\n"
                   "if there is no .pnt file you can try to create an EDF file instead of EDF+.\n\n"
                   "Three files are needed with the extions .eeg, .pnt and .log\n"
                   "this will create an EDF+ file including annotations.\n\n"
                   "If only the .eeg file is available, use: nk2edf -no-annotations <filename>\n"
                   "this will create an EDF file.\n\n",
                   pntfilepath);
            fclose(logfile);
            fclose(inputfile);
            free(log_buf);
            free(sublog_buf);
            return(1);
        }

        rewind(pntfile);
        fread(scratchpad, 16, 1, pntfile);
        scratchpad[16] = 0;
        if(check_device(scratchpad)) {
            printf("error, .pnt file has unknown signature: \"%s\"\n", scratchpad);
            fclose(pntfile);
            fclose(logfile);
            fclose(inputfile);
            free(log_buf);
            free(sublog_buf);
            return(1);
        }
    }

    /*** initialize labels */

    for(i=0; i<256; i++) {
        strcpy(labels[i], "-               ");
    }

    strcpy(labels[0],   "EEG FP1         ");
    strcpy(labels[1],   "EEG FP2         ");
    strcpy(labels[2],   "EEG F3          ");
    strcpy(labels[3],   "EEG F4          ");
    strcpy(labels[4],   "EEG C3          ");
    strcpy(labels[5],   "EEG C4          ");
    strcpy(labels[6],   "EEG P3          ");
    strcpy(labels[7],   "EEG P4          ");
    strcpy(labels[8],   "EEG O1          ");
    strcpy(labels[9],   "EEG O2          ");
    strcpy(labels[10],  "EEG F7          ");
    strcpy(labels[11],  "EEG F8          ");
    strcpy(labels[12],  "EEG T3          ");
    strcpy(labels[13],  "EEG T4          ");
    strcpy(labels[14],  "EEG T5          ");
    strcpy(labels[15],  "EEG T6          ");
    strcpy(labels[16],  "EEG FZ          ");
    strcpy(labels[17],  "EEG CZ          ");
    strcpy(labels[18],  "EEG PZ          ");
    strcpy(labels[19],  "EEG E           ");
    strcpy(labels[20],  "EEG PG1         ");
    strcpy(labels[21],  "EEG PG2         ");
    strcpy(labels[22],  "EEG A1          ");
    strcpy(labels[23],  "EEG A2          ");
    strcpy(labels[24],  "EEG T1          ");
    strcpy(labels[25],  "EEG T2          ");
    for(i=26; i<35; i++) {
        sprintf(labels[i], "EEG X%i          ", i - 25);
    }
    strcpy(labels[35],  "EEG X10         ");
    strcpy(labels[36],  "EEG X11         ");
    for(i=42; i<74; i++) {
        sprintf(labels[i], "DC%02i            ", i - 41);
    }
    strcpy(labels[74],  "EEG BN1         ");
    strcpy(labels[75],  "EEG BN2         ");
    strcpy(labels[76],  "EEG Mark1       ");
    strcpy(labels[77],  "EEG Mark2       ");
    strcpy(labels[100], "EEG X12/BP1     ");
    strcpy(labels[101], "EEG X13/BP2     ");
    strcpy(labels[102], "EEG X14/BP3     ");
    strcpy(labels[103], "EEG X15/BP4     ");
    for(i=104; i<188; i++) {
        sprintf(labels[i], "EEG X%i         ", i - 88);
    }
    for(i=188; i<254; i++) {
        sprintf(labels[i], "EEG X%i        ", i - 88);
    }
    strcpy(labels[255], "Z               ");

    if(read_21e_file(path)) {
        printf("Can not open *.21e file, converter will use default electrode names.\n");
    }

    /***************** start conversion **************************************/

    total_blocks = 0;

    fseeko(inputfile, 0x0091LL, SEEK_SET);
    ctl_block_cnt = fgetc(inputfile);

    printf("ctl_block_cnt = %d\n", ctl_block_cnt);  //***********

    if(ctl_block_cnt==EOF) {
        printf("Error reading inputfile.\n");
        fclose(inputfile);
        if(edfplus) {
            fclose(logfile);
            free(log_buf);
            free(sublog_buf);
        }
        return(1);
    }

    for(i=0; i<ctl_block_cnt; i++) {
        fseeko(inputfile, 0x0092LL + (i * 20LL), SEEK_SET);
        fread((char *)(&ctlblock_address), 4, 1, inputfile);
        fseeko(inputfile, ctlblock_address + 17LL, SEEK_SET);
        datablock_cnt = fgetc(inputfile);

        printf("datablock_cnt = %d\n", datablock_cnt);  //*************

        if(datablock_cnt==EOF) {
            printf("Error reading inputfile.\n");
            fclose(inputfile);
            if(edfplus) {
                fclose(logfile);
                free(log_buf);
                free(sublog_buf);
            }
            return(1);
        }

        for(j=0; j<datablock_cnt; j++) {
            fseeko(inputfile, ctlblock_address + (j * 20LL) + 18LL, SEEK_SET);
            fread((char *)(&wfmblock_address), 4, 1, inputfile);

            printf("wfmblock_address = %d\n", wfmblock_address);  //***********

            /********************************************************************/
            if(edfplus)
                sprintf(path + pathlen - 4, "_%u-%u+.edf", i + 1, j + 1);
            else
                sprintf(path + pathlen - 4, "_%u-%u.edf", i + 1, j + 1);

            outputfile = fopeno(path, "wb");
            if(outputfile==NULL) {
                printf("Can not open file %s for writing.\n", path);
                fclose(inputfile);
                if(edfplus) {
                    fclose(logfile);
                    free(log_buf);
                    free(sublog_buf);
                }
                return(1);
            }

            error = convert_nk2edf(inputfile, outputfile, pntfile, wfmblock_address, edfplus, total_logs, log_buf, read_subevents);
                        
            if(error==1)  printf("Malloc error during conversion\n");
            if(error==2)  printf("Read error during conversion.\n");
            if(error==3)  printf("Write error during conversion.\n");
            if(error==4)  printf("Format error.\n");

            if(fclose(outputfile)) {
                printf("Error closing outputfile.\n");
                fclose(inputfile);
                if(edfplus) {
                    fclose(logfile);
                    fclose(pntfile);
                    free(log_buf);
                    free(sublog_buf);
                }
                return(1);
            }

            if(error) {
                if(edfplus) {
                    fclose(logfile);
                    fclose(pntfile);
                    free(log_buf);
                    free(sublog_buf);
                }
                return(1);
            }

            total_blocks++;

        }
    }

    if(fclose(inputfile))  printf("Error closing inputfile.\n");
    if(edfplus) {
        if(fclose(logfile))  printf("Error closing .log file.\n");
        if(fclose(pntfile))  printf("Error closing .pnt file.\n");
        free(log_buf);
        free(sublog_buf);
    }

    return(0);
}



int convert_nk2edf(FILE *inputfile, FILE *outputfile, FILE *pntfile,  int offset, int edfplus, int n_logs, char *log_buf, int read_subevents)
{
    int i, j, k, p,
        temp,
        channels,
        samplefrequency,
        record_duration,
        raster,
        record_size,
        max_buf_records,
        bufsize,
        records_in_buf,
        seconds,
        deci_seconds,
        left_records,
        elapsed_time,
        error,
        n_log_processed;

    char *buf,
         *annotations,
         scratchpad[48];

    /* filter events */

    printf("n_logs = %d\n", n_logs);  //*************

    for(i=0; i<n_logs; i++) {
        elapsed_time = 36000 * (log_buf[(i * 45) + 20] - 48);
        elapsed_time += 3600 * (log_buf[(i * 45) + 21] - 48);
        elapsed_time += 600 * (log_buf[(i * 45) + 22] - 48);
        elapsed_time += 60 * (log_buf[(i * 45) + 23] - 48);
        elapsed_time += 10 * (log_buf[(i * 45) + 24] - 48);
        elapsed_time += log_buf[(i * 45) + 25] - 48;
        if(elapsed_time>=total_elapsed_time) break;
    }

    printf("elapsed_time = %d, total_elapsed_time = %d\n", elapsed_time, total_elapsed_time);

    log_buf += i * 45;
    n_logs -= i;

    /************************* write EDF-header ***************************************/

    rewind(outputfile);

    fprintf(outputfile, "0       ");

    if(edfplus) {
        error = 0;
        fseeko(pntfile, 0x0604LL, SEEK_SET);
        fread(scratchpad, 10, 1, pntfile);
        scratchpad[10] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<10; i++) {
            if(scratchpad[i]==0)  break;
            if(scratchpad[i]==' ')  scratchpad[i] = '_';
        }
        if(i) {
            p = i;
            fwrite(scratchpad, i, 1, outputfile);
        } else {
            fputc('X', outputfile);
            p = 1;
        }
        fputc(' ', outputfile);
        p++;

        fseeko(pntfile, 0x064aLL, SEEK_SET);
        fread(scratchpad, 6, 1, pntfile);
        if(!strncmp(scratchpad, "Male", 4))  fputc('M', outputfile);
        else {
            if(!strncmp(scratchpad, "Female", 6))  fputc('F', outputfile);
            else  fputc('X', outputfile);
        }
        p++;
        fputc(' ', outputfile);
        p++;

        fseeko(pntfile, 0x0668LL, SEEK_SET);
        fread(scratchpad, 2, 1, pntfile);
        scratchpad[2] = 0;
        temp = atoi(scratchpad);
        if((temp<1)||(temp>31))  error = 1;
        for(i=0; i<2; i++) {
            if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                error = 1;
                break;
            }
        }
        fseeko(pntfile, 0x0665LL, SEEK_SET);
        fread(scratchpad, 2, 1, pntfile);
        scratchpad[2] = 0;
        temp = atoi(scratchpad);
        if((temp<1)||(temp>12))  error = 1;
        fseeko(pntfile, 0x0660LL, SEEK_SET);
        fread(scratchpad, 4, 1, pntfile);
        scratchpad[4] = 0;
        temp = atoi(scratchpad);
        if((temp<1)||(temp>9999))  error = 1;
        for(i=0; i<4; i++) {
            if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                error = 1;
                break;
            }
        }

        if(error) {
            fputc('X', outputfile);
            p++;
        } else {
            fseeko(pntfile, 0x0668LL, SEEK_SET);
            fread(scratchpad, 2, 1, pntfile);
            scratchpad[2] = 0;
            temp = atoi(scratchpad);
            if((temp<1)||(temp>31)) {
                sprintf(scratchpad, "01");
                error = 1;
            }
            for(i=0; i<2; i++) {
                if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                    sprintf(scratchpad, "01");
                    error = 1;
                    break;
                }
            }
            fwrite(scratchpad, 2, 1, outputfile);
            p += 2;
            fputc('-', outputfile);
            p++;
            fseeko(pntfile, 0x0665LL, SEEK_SET);
            fread(scratchpad, 2, 1, pntfile);
            scratchpad[2] = 0;
            temp = atoi(scratchpad);
            switch(temp) {
            case  1 :
                fwrite("JAN", 3, 1, outputfile);
                break;
            case  2 :
                fwrite("FEB", 3, 1, outputfile);
                break;
            case  3 :
                fwrite("MAR", 3, 1, outputfile);
                break;
            case  4 :
                fwrite("APR", 3, 1, outputfile);
                break;
            case  5 :
                fwrite("MAY", 3, 1, outputfile);
                break;
            case  6 :
                fwrite("JUN", 3, 1, outputfile);
                break;
            case  7 :
                fwrite("JUL", 3, 1, outputfile);
                break;
            case  8 :
                fwrite("AUG", 3, 1, outputfile);
                break;
            case  9 :
                fwrite("SEP", 3, 1, outputfile);
                break;
            case 10 :
                fwrite("OCT", 3, 1, outputfile);
                break;
            case 11 :
                fwrite("NOV", 3, 1, outputfile);
                break;
            case 12 :
                fwrite("DEC", 3, 1, outputfile);
                break;
            default :
                fwrite("JAN", 3, 1, outputfile);
                error = 1;
                break;
            }
            p += 3;
            fputc('-', outputfile);
            p++;
            fseeko(pntfile, 0x0660LL, SEEK_SET);
            fread(scratchpad, 4, 1, pntfile);
            scratchpad[4] = 0;
            temp = atoi(scratchpad);
            if((temp<1)||(temp>9999)) {
                sprintf(scratchpad, "1800");
                error = 1;
            }
            for(i=0; i<4; i++) {
                if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                    sprintf(scratchpad, "1800");
                    error = 1;
                    break;
                }
            }
            fwrite(scratchpad, 4, 1, outputfile);
            p += 4;
        }

        fputc(' ', outputfile);
        p++;

        fseeko(pntfile, 0x062eLL, SEEK_SET);
        fread(scratchpad, 20, 1, pntfile);
        scratchpad[20] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<20; i++) {
            if(scratchpad[i]==0)  break;
            if(scratchpad[i]==' ')  scratchpad[i] = '_';
        }
        if(i) {
            p += i;
            fwrite(scratchpad, i, 1, outputfile);
        } else {
            fputc('X', outputfile);
            p++;
        }

        for(i=0; i<80-p; i++)  fputc(' ', outputfile);

        fwrite("Startdate ", 10, 1, outputfile);
        p = 10;
        error = 0;
        fseeko(pntfile, 0x0046LL, SEEK_SET);
        fread(scratchpad, 2, 1, pntfile);
        scratchpad[2] = 0;
        temp = atoi(scratchpad);
        if((temp<1)||(temp>31))  error = 1;
        for(i=0; i<2; i++) {
            if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                error = 1;
                break;
            }
        }
        fseeko(pntfile, 0x0044LL, SEEK_SET);
        fread(scratchpad, 2, 1, pntfile);
        scratchpad[2] = 0;
        temp = atoi(scratchpad);
        if((temp<1)||(temp>12))  error = 1;
        fseeko(pntfile, 0x0040LL, SEEK_SET);
        fread(scratchpad, 4, 1, pntfile);
        scratchpad[4] = 0;
        temp = atoi(scratchpad);
        if((temp<1)||(temp>9999))  error = 1;
        for(i=0; i<4; i++) {
            if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                error = 1;
                break;
            }
        }

        if(error) {
            fputc('X', outputfile);
            p++;
        } else {
            fseeko(pntfile, 0x0046LL, SEEK_SET);
            fread(scratchpad, 2, 1, pntfile);
            scratchpad[2] = 0;
            temp = atoi(scratchpad);
            if((temp<1)||(temp>31))  sprintf(scratchpad, "01");
            for(i=0; i<2; i++) {
                if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                    sprintf(scratchpad, "01");
                    break;
                }
            }
            fwrite(scratchpad, 2, 1, outputfile);
            fputc('-', outputfile);
            fseeko(pntfile, 0x0044LL, SEEK_SET);
            fread(scratchpad, 2, 1, pntfile);
            scratchpad[2] = 0;
            temp = atoi(scratchpad);
            switch(temp) {
            case  1 :
                fwrite("JAN", 3, 1, outputfile);
                break;
            case  2 :
                fwrite("FEB", 3, 1, outputfile);
                break;
            case  3 :
                fwrite("MAR", 3, 1, outputfile);
                break;
            case  4 :
                fwrite("APR", 3, 1, outputfile);
                break;
            case  5 :
                fwrite("MAY", 3, 1, outputfile);
                break;
            case  6 :
                fwrite("JUN", 3, 1, outputfile);
                break;
            case  7 :
                fwrite("JUL", 3, 1, outputfile);
                break;
            case  8 :
                fwrite("AUG", 3, 1, outputfile);
                break;
            case  9 :
                fwrite("SEP", 3, 1, outputfile);
                break;
            case 10 :
                fwrite("OCT", 3, 1, outputfile);
                break;
            case 11 :
                fwrite("NOV", 3, 1, outputfile);
                break;
            case 12 :
                fwrite("DEC", 3, 1, outputfile);
                break;
            default :
                fwrite("JAN", 3, 1, outputfile);
                break;
            }
            fputc('-', outputfile);
            fseeko(pntfile, 0x0040LL, SEEK_SET);
            fread(scratchpad, 4, 1, pntfile);
            scratchpad[4] = 0;
            temp = atoi(scratchpad);
            if((temp<1)||(temp>9999))  sprintf(scratchpad, "1800");
            for(i=0; i<4; i++) {
                if((scratchpad[i]<'0')||(scratchpad[i]>'9')) {
                    sprintf(scratchpad, "1800");
                    break;
                }
            }
            fwrite(scratchpad, 4, 1, outputfile);
            p += 11;
        }

        fputc(' ', outputfile);
        p++;

        fseeko(pntfile, 0x061cLL, SEEK_SET);
        fread(scratchpad, 10, 1, pntfile);
        scratchpad[10] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<10; i++) {
            if(scratchpad[i]==0)  break;
            if(scratchpad[i]==' ')  scratchpad[i] = '_';
        }
        if(i) {
            p += i;
            fwrite(scratchpad, i, 1, outputfile);
        } else {
            fputc('X', outputfile);
            p++;
        }

        fputc(' ', outputfile);
        p++;

        fseeko(pntfile, 0x06aaLL, SEEK_SET);
        fread(scratchpad, 20, 1, pntfile);
        scratchpad[20] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<20; i++) {
            if(scratchpad[i]==0)  break;
            if(scratchpad[i]==' ')  scratchpad[i] = '_';
        }
        if(i) {
            p += i;
            fwrite(scratchpad, i, 1, outputfile);
        } else {
            fputc('X', outputfile);
            p++;
        }

        fputc(' ', outputfile);
        p++;

        fwrite("Nihon_Kohden_", 13, 1, outputfile);
        p += 13;
        rewind(inputfile);
        fread(scratchpad, 16, 1, inputfile);
        scratchpad[16] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<16; i++) {
            if(scratchpad[i]==0)  break;
            if(scratchpad[i]==' ')  scratchpad[i] = '_';
        }
        fwrite(scratchpad, i, 1, outputfile);
        p += i;

        for(i=0; i<80-p; i++)  fputc(' ', outputfile);
    } else {
        fseeko(inputfile, 0x004fLL, SEEK_SET);
        fread(scratchpad, 32, 1, inputfile);
        scratchpad[32] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<32; i++) {
            if(scratchpad[i]==0)  break;
        }
        p = 80 - i;
        fwrite(scratchpad, i, 1, outputfile);
        for(i=0; i<p; i++)  fputc(' ', outputfile);

        fwrite("Nihon Kohden ", 13, 1, outputfile);
        rewind(inputfile);
        fread(scratchpad, 16, 1, inputfile);
        scratchpad[16] = 0;
        latin1_to_ascii(scratchpad, strlen(scratchpad));
        for(i=0; i<16; i++) {
            if(scratchpad[i]==0)  break;
        }
        p = 67 - i;
        fwrite(scratchpad, i, 1, outputfile);
        for(i=0; i<p; i++)  fputc(' ', outputfile);
    }

    fseeko(inputfile, 0x0016LL + offset, SEEK_SET);
    temp = fgetc(inputfile);
    fprintf(outputfile, "%02u.", ((temp >> 4) * 10) + (temp & 15));
    fseeko(inputfile, 0x0015LL + offset, SEEK_SET);
    temp = fgetc(inputfile);
    fprintf(outputfile, "%02u.", ((temp >> 4) * 10) + (temp & 15));
    fseeko(inputfile, 0x0014LL + offset, SEEK_SET);
    temp = fgetc(inputfile);
    fprintf(outputfile, "%02u", ((temp >> 4) * 10) + (temp & 15));

    fseeko(inputfile, 0x0017LL + offset, SEEK_SET);
    temp = fgetc(inputfile);
    fprintf(outputfile, "%02u.", ((temp >> 4) * 10) + (temp & 15));
    temp = fgetc(inputfile);
    fprintf(outputfile, "%02u.", ((temp >> 4) * 10) + (temp & 15));
    temp = fgetc(inputfile);
    fprintf(outputfile, "%02u", ((temp >> 4) * 10) + (temp & 15));

    fseeko(inputfile, 0x0026LL + offset, SEEK_SET);
    channels = fgetc(inputfile) + 1;
    if(edfplus) {
        fprintf(outputfile, "%-8u", (channels + 1) * 256 + 256);
        fprintf(outputfile, "EDF+C");
        for(i=0; i<39; i++)  fputc(' ', outputfile);
    } else {
        fprintf(outputfile, "%-8u", channels * 256 + 256);
        for(i=0; i<44; i++)  fputc(' ', outputfile);
    }
    fseeko(inputfile, 0x001cLL + offset, SEEK_SET);
    fread((char *)(&record_duration), 4, 1, inputfile);
    if((record_duration < 10) || (record_duration > 99999999)) {
        return(4);
    }
    fprintf(outputfile, "%-8u", record_duration);
    fprintf(outputfile, "0.1     ");
    if(edfplus)  fprintf(outputfile, "%-4u", channels + 1);
    else  fprintf(outputfile, "%-4u", channels);

    for(i=0; i<(channels - 1); i++) {
        fseeko(inputfile, (long long)(0x0027 + (i * 10) + offset), SEEK_SET);
        temp = fgetc(inputfile);
        if((temp < 0) || (temp > 255)) {
            fprintf(outputfile, "-               ");
        } else {
            fprintf(outputfile, "%s", labels[temp]);
        }
    }

    fprintf(outputfile, "Events/Markers  ");

    if(edfplus)  fprintf(outputfile, "EDF Annotations ");

    for(i=0; i<(channels * 80); i++)  fputc(' ', outputfile);
    if(edfplus)  for(i=0; i<80; i++)  fputc(' ', outputfile);

    for(i=0; i<(channels - 1); i++) {
        fseeko(inputfile, 0x0027LL + (i * 10LL) + offset, SEEK_SET);
        temp = fgetc(inputfile);
        if(((temp<42)||(temp>73)) && (temp!=76) && (temp!=77))  fprintf(outputfile, "uV      ");
        else  fprintf(outputfile, "mV      ");
    }
    fprintf(outputfile, "        ");
    if(edfplus)  fprintf(outputfile, "        ");

    for(i=0; i<(channels - 1); i++) {
        fseeko(inputfile, 0x0027LL + (i * 10LL) + offset, SEEK_SET);
        temp = fgetc(inputfile);
        if(((temp<42)||(temp>73)) && (temp!=76) && (temp!=77))  fprintf(outputfile, "-3200   ");
        else  fprintf(outputfile, "-12002.9");
    }
    fprintf(outputfile, "-1      ");
    if(edfplus)  fprintf(outputfile, "-1      ");

    for(i=0; i<(channels - 1); i++) {
        fseeko(inputfile, 0x0027LL + (i * 10LL) + offset, SEEK_SET);
        temp = fgetc(inputfile);
        if(((temp<42)||(temp>73)) && (temp!=76) && (temp!=77))  fprintf(outputfile, "3199.902");
        else  fprintf(outputfile, "12002.56");
    }
    fprintf(outputfile, "1       ");
    if(edfplus)  fprintf(outputfile, "1       ");

    for(i=0; i<channels; i++)  fprintf(outputfile, "-32768  ");
    if(edfplus)  fprintf(outputfile, "-32768  ");

    for(i=0; i<channels; i++)  fprintf(outputfile, "32767   ");
    if(edfplus)  fprintf(outputfile, "32767   ");

    for(i=0; i<(channels * 80); i++)  fputc(' ', outputfile);
    if(edfplus)  for(i=0; i<80; i++)  fputc(' ', outputfile);

    fseeko(inputfile, 0x001bLL + offset, SEEK_SET);
    samplefrequency = fgetc(inputfile) * 256;
    fseeko(inputfile, 0x001aLL + offset, SEEK_SET);
    samplefrequency += fgetc(inputfile);
    samplefrequency &= 0x3fff;
    for(i=0; i<channels; i++)  fprintf(outputfile, "%-8u", samplefrequency / 10);
    if(edfplus)  fprintf(outputfile, "%-8u", ANNOT_TRACKSIZE / 2);

    for(i=0; i<(channels * 32); i++)  fputc(' ', outputfile);
    if(edfplus)  for(i=0; i<32; i++)  fputc(' ', outputfile);

    /************************* write data ****************************************************/

    bufsize = 4194304;
    buf = (char *)malloc(bufsize);
    if(buf==NULL)  
        return(1);

    record_size = (samplefrequency / 10) * channels * 2;

    printf("samplefrequency = %d, channels = %d\n", samplefrequency, channels);  //*****
    printf("record_size = %d\n", record_size);  //***********

    if(edfplus)  
        record_size += ANNOT_TRACKSIZE;
    
    printf("record_size = %d\n", record_size);  //***********

    max_buf_records = bufsize / record_size;

    raster = (samplefrequency / 10) * 2;

    printf("max_buf_records = %d, raster = %d\n", max_buf_records, raster);  //*******

    seconds = 0;
    deci_seconds = 0;
    n_log_processed = 0;

    fseeko(inputfile, (long long)(0x0027LL + offset + ((channels - 1) * 10LL)), SEEK_SET);

    left_records = record_duration;

    printf("record_duration = %d\n", record_duration);  //**********

    while(left_records) {
        if(left_records>max_buf_records)  records_in_buf = max_buf_records;
        else  records_in_buf = left_records;

        for(i=0; i<records_in_buf; i++) {
            for(j=0; j<raster; j+=2) {
                for(k=0; k<(channels - 1); k++) {
                    buf[j+(k*raster)+(i*record_size)] = fgetc(inputfile);
                    buf[j+(k*raster)+(i*record_size)+1] = fgetc(inputfile) + 128;

                    // putchar(buf[j+(k*raster)+(i*record_size)]);  //************

                }
                buf[j+(k*raster)+(i*record_size)] = fgetc(inputfile);
                temp = fgetc(inputfile);
                if(temp==EOF) {
                    free(buf);
                    return(2);
                }
                buf[j+(k*raster)+(i*record_size)+1] = temp;
            }
            if(edfplus) {
                annotations = buf + (i * record_size) + (raster * channels);
                memset(annotations, 0, ANNOT_TRACKSIZE);
                p = sprintf(annotations, "%+i.%i", seconds, deci_seconds);
                annotations[p++] = 20;
                annotations[p++] = 20;
                for( ; n_log_processed < n_logs; n_log_processed++) {
                    elapsed_time = 36000 * (log_buf[(n_log_processed * 45) + 20] - 48);
                    elapsed_time += 3600 * (log_buf[(n_log_processed * 45) + 21] - 48);
                    elapsed_time += 600 * (log_buf[(n_log_processed * 45) + 22] - 48);
                    elapsed_time += 60 * (log_buf[(n_log_processed * 45) + 23] - 48);
                    elapsed_time += 10 * (log_buf[(n_log_processed * 45) + 24] - 48);
                    elapsed_time += log_buf[(n_log_processed * 45) + 25] - 48;
                    if(elapsed_time>=total_elapsed_time) {
                        elapsed_time -= total_elapsed_time;
                        if(elapsed_time<(record_duration / 10)) {
                            p++;
                            p += sprintf(annotations + p, "%+i", elapsed_time);
                            if(read_subevents) {
                                annotations[p] = '.';
                                p++;
                                strncpy(annotations + p, log_buf + (n_log_processed * 45) + 26, 3);
                                p += 3;
                            }
                            annotations[p++] = 20;
                            strncpy(annotations + p, log_buf + (n_log_processed * 45), 20);
                            p += 20;
                            annotations[p] = 20;

                            n_log_processed++;

                            break;
                        }
                    }
                }
            }
            if(++deci_seconds>9) {
                deci_seconds = 0;
                seconds++;
            }
        }

        

        if(fwrite(buf, records_in_buf * record_size, 1, outputfile)!=1) {
            free(buf);
            return(3);
        }
        left_records -= records_in_buf;
    }

    total_elapsed_time += record_duration / 10;

    free(buf);

    return(0);
}


int check_device(char *str)
{
    int error = 1;

    if(!strncmp(str, "EEG-1100A V01.00", 16))  error = 0;
    if(!strncmp(str, "EEG-1100B V01.00", 16))  error = 0;
    if(!strncmp(str, "EEG-1100C V01.00", 16))  error = 0;
    if(!strncmp(str, "QI-403A   V01.00", 16))  error = 0;
    if(!strncmp(str, "QI-403A   V02.00", 16))  error = 0;
    if(!strncmp(str, "EEG-2100  V01.00", 16))  error = 0;
    if(!strncmp(str, "EEG-2100  V02.00", 16))  error = 0;
    if(!strncmp(str, "DAE-2100D V01.30", 16))  error = 0;
    if(!strncmp(str, "DAE-2100D V02.00", 16))  error = 0;
    if(!strncmp(str, "EEG-1100A V02.00", 16))  error = 0;
    if(!strncmp(str, "EEG-1100B V02.00", 16))  error = 0;
    if(!strncmp(str, "EEG-1100C V02.00", 16))  error = 0;
    /* workaround for log file quirk where the last character of the version string is missing: */
    if((!strncmp(str, "EEG-1100A V02.0",  15)) && (str[15] == 0))  error = 0;

    return(error);
}




void latin1_to_utf8(char *latin1_str, int len)
{
    int i, j;

    unsigned char *str, tmp_str[1024];


    str = (unsigned char *)latin1_str;

    j = 0;

    for(i=0; i<len; i++) {
        tmp_str[j] = str[i];

        if(str[i]==0)  break;

        if(str[i]<32) tmp_str[j] = '.';

        if((str[i]>126)&&(str[i]<160))  tmp_str[j] = '.';

        if(str[i]>159) {
            if((len-j)<2) {
                tmp_str[j] = ' ';
            } else {
                tmp_str[j] = 192 + (str[i]>>6);
                j++;
                tmp_str[j] = 128 + (str[i]&63);
            }
        }

        j++;

        if(j>=len)  break;
    }

    for(i=0; i<len; i++) {
        str[i] = tmp_str[i];

        if(str[i]==0)  return;
    }
}



void latin1_to_ascii(char *str, int len)
{
    int i, value;

    for(i=0; i<len; i++) {
        value = *((unsigned char *)(str + i));

        if((value>31)&&(value<127)) {
            continue;
        }

        switch(value) {
        case 128 :
            str[i] = 'E';
            break;

        case 130 :
            str[i] = ',';
            break;

        case 131 :
            str[i] = 'F';
            break;

        case 132 :
            str[i] = '\"';
            break;

        case 133 :
            str[i] = '.';
            break;

        case 134 :
            str[i] = '+';
            break;

        case 135 :
            str[i] = '+';
            break;

        case 136 :
            str[i] = '^';
            break;

        case 137 :
            str[i] = 'm';
            break;

        case 138 :
            str[i] = 'S';
            break;

        case 139 :
            str[i] = '<';
            break;

        case 140 :
            str[i] = 'E';
            break;

        case 142 :
            str[i] = 'Z';
            break;

        case 145 :
            str[i] = '`';
            break;

        case 146 :
            str[i] = '\'';
            break;

        case 147 :
            str[i] = '\"';
            break;

        case 148 :
            str[i] = '\"';
            break;

        case 149 :
            str[i] = '.';
            break;

        case 150 :
            str[i] = '-';
            break;

        case 151 :
            str[i] = '-';
            break;

        case 152 :
            str[i] = '~';
            break;

        case 154 :
            str[i] = 's';
            break;

        case 155 :
            str[i] = '>';
            break;

        case 156 :
            str[i] = 'e';
            break;

        case 158 :
            str[i] = 'z';
            break;

        case 159 :
            str[i] = 'Y';
            break;

        case 171 :
            str[i] = '<';
            break;

        case 180 :
            str[i] = '\'';
            break;

        case 181 :
            str[i] = 'u';
            break;

        case 187 :
            str[i] = '>';
            break;

        case 191 :
            str[i] = '\?';
            break;

        case 192 :
            str[i] = 'A';
            break;

        case 193 :
            str[i] = 'A';
            break;

        case 194 :
            str[i] = 'A';
            break;

        case 195 :
            str[i] = 'A';
            break;

        case 196 :
            str[i] = 'A';
            break;

        case 197 :
            str[i] = 'A';
            break;

        case 198 :
            str[i] = 'E';
            break;

        case 199 :
            str[i] = 'C';
            break;

        case 200 :
            str[i] = 'E';
            break;

        case 201 :
            str[i] = 'E';
            break;

        case 202 :
            str[i] = 'E';
            break;

        case 203 :
            str[i] = 'E';
            break;

        case 204 :
            str[i] = 'I';
            break;

        case 205 :
            str[i] = 'I';
            break;

        case 206 :
            str[i] = 'I';
            break;

        case 207 :
            str[i] = 'I';
            break;

        case 208 :
            str[i] = 'D';
            break;

        case 209 :
            str[i] = 'N';
            break;

        case 210 :
            str[i] = 'O';
            break;

        case 211 :
            str[i] = 'O';
            break;

        case 212 :
            str[i] = 'O';
            break;

        case 213 :
            str[i] = 'O';
            break;

        case 214 :
            str[i] = 'O';
            break;

        case 215 :
            str[i] = 'x';
            break;

        case 216 :
            str[i] = 'O';
            break;

        case 217 :
            str[i] = 'U';
            break;

        case 218 :
            str[i] = 'U';
            break;

        case 219 :
            str[i] = 'U';
            break;

        case 220 :
            str[i] = 'U';
            break;

        case 221 :
            str[i] = 'Y';
            break;

        case 222 :
            str[i] = 'I';
            break;

        case 223 :
            str[i] = 's';
            break;

        case 224 :
            str[i] = 'a';
            break;

        case 225 :
            str[i] = 'a';
            break;

        case 226 :
            str[i] = 'a';
            break;

        case 227 :
            str[i] = 'a';
            break;

        case 228 :
            str[i] = 'a';
            break;

        case 229 :
            str[i] = 'a';
            break;

        case 230 :
            str[i] = 'e';
            break;

        case 231 :
            str[i] = 'c';
            break;

        case 232 :
            str[i] = 'e';
            break;

        case 233 :
            str[i] = 'e';
            break;

        case 234 :
            str[i] = 'e';
            break;

        case 235 :
            str[i] = 'e';
            break;

        case 236 :
            str[i] = 'i';
            break;

        case 237 :
            str[i] = 'i';
            break;

        case 238 :
            str[i] = 'i';
            break;

        case 239 :
            str[i] = 'i';
            break;

        case 240 :
            str[i] = 'd';
            break;

        case 241 :
            str[i] = 'n';
            break;

        case 242 :
            str[i] = 'o';
            break;

        case 243 :
            str[i] = 'o';
            break;

        case 244 :
            str[i] = 'o';
            break;

        case 245 :
            str[i] = 'o';
            break;

        case 246 :
            str[i] = 'o';
            break;

        case 247 :
            str[i] = '-';
            break;

        case 248 :
            str[i] = '0';
            break;

        case 249 :
            str[i] = 'u';
            break;

        case 250 :
            str[i] = 'u';
            break;

        case 251 :
            str[i] = 'u';
            break;

        case 252 :
            str[i] = 'u';
            break;

        case 253 :
            str[i] = 'y';
            break;

        case 254 :
            str[i] = 't';
            break;

        case 255 :
            str[i] = 'y';
            break;

        default  :
            str[i] = ' ';
            break;
        }
    }
}


int read_21e_file(char *e21filepath)
{
    int n,
        flag_eleclines=0,
        idx;

    char *electrode_name,
         electrode_name_buffer[ELECTRODE_NAME_MAXLEN],
         scratchpad[64],
         *charpntr;

    FILE *inputfile;


    e21filepath[strlen(e21filepath) - 4] = 0;
    strcat(e21filepath, ".21E");
    inputfile = fopeno(e21filepath, "rb");
    if(inputfile==NULL) {
        e21filepath[strlen(e21filepath) - 4] = 0;
        strcat(e21filepath, ".21e");
        inputfile = fopeno(e21filepath, "rb");
        if(inputfile==NULL) {
            return(1);
        }
    }

    while (!feof(inputfile)) {
        charpntr = fgets(electrode_name_buffer, ELECTRODE_NAME_MAXLEN-1, inputfile);

        if(charpntr == NULL) {
            break;
        }

        if(strncmp(electrode_name_buffer, ELECTRODE_TAG, strlen(ELECTRODE_TAG)) == 0) {
            flag_eleclines = 1;
        } else {
            if(strncmp(electrode_name_buffer, ELECTRODE_UNTAG, strlen(ELECTRODE_UNTAG)) == 0) {
                flag_eleclines = 0;
            }
        }

        if(flag_eleclines) {
            if(strtok(electrode_name_buffer, "=") != NULL) {
                idx = atoi(electrode_name_buffer);

                electrode_name = strtok(NULL, "=");

                if(electrode_name != NULL) {
                    n = strlen(electrode_name);

                    if((n > 0)&&(electrode_name[n-1] == 10)) {
                        electrode_name[n-1] = 0;
                    }

                    if((n > 1)&&(electrode_name[n-2] == 13)) {
                        electrode_name[n-2] = 0;
                    }

                    n = strlen(electrode_name);

                    if((idx >= 0) && (idx < 256)) {
                        if(n > 0) {
                            strncpy(scratchpad, electrode_name, 16);

                            strcat(scratchpad, "                ");

                            latin1_to_ascii(scratchpad, 16);

                            scratchpad[16] = 0;

                            strcpy(labels[idx], scratchpad);
                        } else {
                            strcpy(labels[idx], "-               ");
                        }
                    }
                }
            }
        }
    }

    fclose(inputfile);

    return(0);
}


