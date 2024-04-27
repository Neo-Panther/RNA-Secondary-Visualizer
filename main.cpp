#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;
/**Tests
 * ACCGGUAGU = 2
 * CCCCCGAAAAAU = 2
 * 
 * RNA Central Test Sequences:
 * GGCGAAGCCCGCCUGUGCGGGCUA = 7 - Perfect Match
 * CUUGCUGAGGUGCACACAGCAAG = ours 8, theirs 7: they miss pair between 11(U) and 16(A)
 * AGAUCUGAGCCUGGGAGCUCUCU = 7 - Perfect Match
 * AAACCGUUACCAUUACUGAGUUU = ours 8, theirs 1: they miss multiple possible pairings 
 * GGCGAAGCCCGCCUGUGCGGGCUA = 7 - Perfect Match
 * 
*/

int main(int argc, char* argv[]){
  if(argc < 2){
    printf("Please enter rna sequence after the file name.\n");
    return 0;
  }
  setbuf(stdout, NULL);
  char *rna = argv[1]; int n, m;
  n = strlen(rna);
  int i = 0;
  for(i = 0; rna[i] == ' '; i++);
  for(int j = 0; j < n; j++) rna[j] = rna[i+j];
  for(i = n-2; i >= 0 && rna[i] == ' '; i--);
  rna[i+2] = 0;
  n = strlen(rna);
  char* out = (char*)calloc(n+1, sizeof(char));
  for(int i = 0; i < n; i++){
    out[i] = '.';
  }
  out[n] = 0;
  if(n > 5){
    int** to_clear = (int**)calloc(n*n, sizeof(int*));
    int to_clear_i = 0;
    int**** OPTa = (int****)calloc(n-5, sizeof(int***));
    int** OPT = (int**)calloc(n-5, sizeof(int*));
    for(int i = 0; i < n-5; i++){
      OPT[i] = (int*)calloc(n-5, sizeof(int));
      OPTa[i] = (int***)calloc(n-5, sizeof(int**));
      for(int j = 0; j < n-5; j++){
        OPTa[i][j] = (int**)calloc(n, sizeof(int*));
      }
    }
    for(int k = 5; k < n; k++){
      for(int i = 0; i < n-k; i++){
        int j = i+k-5;
        if(j > 0){
          OPT[i][j] = OPT[i][j-1];
          for(int l = 0; l < OPT[i][j]; l++){
            OPTa[i][j][l] = OPTa[i][j-1][l];
          }
        } else {
          OPT[i][j] = 0;
        }
        for(int t = i; t <= j; t++){
          if((rna[t] == 'A' && rna[j+5] == 'U') || (rna[t] == 'C' && rna[j+5] == 'G') || (rna[t] == 'U' && rna[j+5] == 'A') || (rna[t] == 'G' && rna[j+5] == 'C')){
            int topt = 1;
            if(t-6 >= 0) topt += OPT[i][t-6];
            if(t < n-6 && j > 0) topt += OPT[t+1][j-1];
            if(OPT[i][j] <= topt){
              int* tmp = (int*)calloc(2, sizeof(int));
              tmp[0] = t+1;
              tmp[1] = j+6;
              to_clear[to_clear_i++] = tmp;
              int l1 = 0;
              OPTa[i][j][l1++] = tmp;
              if(t-6 >= 0){
                for(int l = 0; l < OPT[i][t-6]; l++){
                  OPTa[i][j][l1++] = OPTa[i][t-6][l];
                }
              }
              if(t < n-6 && j > 0) {
                for(int l = 0; l < OPT[t+1][j-1]; l++){
                  OPTa[i][j][l1++] = OPTa[t+1][j-1][l];
                }
              }
              OPT[i][j] = topt;
            }
          }
        }
      }
    }
    printf("Answer: %d Pairs Matched\n", OPT[0][n-6]);
    m = OPT[0][n-6];
    printf("::Pairs matched::\n");
    for(int i=0;i<m;i++){
      cout<<i+1<<". "<<OPTa[0][n-6][i][0]<<" "<<OPTa[0][n-6][i][1]<<endl;
      if(OPTa[0][n-6][i][0] > OPTa[0][n-6][i][1]){
        out[OPTa[0][n-6][i][0]-1] = ')';
        out[OPTa[0][n-6][i][1]-1] = '(';
      } else {
        out[OPTa[0][n-6][i][0]-1] = '(';
        out[OPTa[0][n-6][i][1]-1] = ')';
      }
    }
    // free created arrays
    for(int i = 0; i < n-5; i++){
      for(int j = 0; j < n-5; j++){
        free(OPTa[i][j]);
      }
      free(OPT[i]);
      free(OPTa[i] = (int***)calloc(n-5, sizeof(int**)));
    }
    free(OPT);
    free(OPTa);
  } else {
    m = 0;
    printf("Answer: 0 Pairs Matched\n");
    printf("No Pairs Matched\n");
  }
  // compile file to svg
  char command[1000];
  sprintf(command, "java -cp assets/bin/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN %s -o viz.svg -algorithm naview", rna, out);
  int result = system(command);
  printf("Vizualization saved in 'viz.svg'\n");
  free(out);
  return 0;
}