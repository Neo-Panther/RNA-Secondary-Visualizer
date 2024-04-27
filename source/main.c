/**
 * @file main.c
 * @author Aryan Gupta
 * @brief The file containing the implementation of the RNA secondary structure algorithm
 * @version 1.0
 * @date 2024-04-27
 *  
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
// Access the Windows API (for system calls)
#include <windows.h>

/* Tests
 * ACCGGUAGU = 2
 * CCCCCGAAAAAU = 2
 * UUUUUUUUUUUU = 0
 * RNA Central Test Sequences:
 * GGCGAAGCCCGCCUGUGCGGGCUA = 7 - Perfect Match (URS00006569A0)
 * CUUGCUGAGGUGCACACAGCAAG = ours 8, theirs 7: they miss pair between 11(U) and 16(A) (URS000080E21C)
 * AGAUCUGAGCCUGGGAGCUCUCU = 7 - Perfect Match (URS0000A76371)
 * AAACCGUUACCAUUACUGAGUUU = ours 8, theirs 1: they miss multiple possible pairings (URS00002AB26F)
 * CUUUCAAUCCUCUUCUUGAGAUUC = 7 - Different pairs: same bases (URS000064F01F)
 * 
*/

/**
 * @brief Trims white spaces from start and end of the input string.
 * 
 * @param [in out] str The input string
 */
void lrtrim(char* str){
  int n = strlen(str);
  int i = 0;
  for(i = 0; str[i] == ' '; i++);
  for(int j = 0; j < n; j++) str[j] = str[i+j];
  for(i = n-2; i >= 0 && str[i] == ' '; i--);
  str[i+2] = 0;
}

/**
 * @brief The Function containg the algorithm implementation.
 * 
 * @param [out] to_clear Stores pointers of base pairs to be freed
 * @param [out] to_clear_i Index of the next entry in to_clear array
 * @param [out] OPTa Dynamic Programming Array, stores the pairs maximizing the maximum number of pairs found till now
 * @param [out] OPT Dynamic Programming Array, stores maximum number of pairs found till now
 * @param [in] n Length of the input RNA sequence
 * @param [in] rna The rna sequence
 * @return int 
 */
int rnaFold(int** to_clear, int to_clear_i, int**** OPTa, int** OPT, int n, char* rna){
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
  return to_clear_i;
}

/**
 * @brief Builds the visualization output from the algorithm's output
 * 
 * @param [in] m Number of pairs in algorithms output
 * @param [in] OPTa Output pairs from the algorithm
 * @param [out] out Output to be used in visualization
 * @param [in] n Number of bases in input RNA
 */
void buildOut(int m, int**** OPTa, char* out, int n){
  for(int i=0;i<m;i++){
    printf("%d. %d - %d\n", i+1, OPTa[0][n-6][i][0], OPTa[0][n-6][i][1]);
    if(OPTa[0][n-6][i][0] > OPTa[0][n-6][i][1]){
      out[OPTa[0][n-6][i][0]-1] = ')';
      out[OPTa[0][n-6][i][1]-1] = '(';
    } else {
      out[OPTa[0][n-6][i][0]-1] = '(';
      out[OPTa[0][n-6][i][1]-1] = ')';
    }
  }
}

/**
 * @brief Main Entry point of the program.
 * 
 * @param [in] argc number of input arguments
 * @param [in] argv the input arguments
 * @return int exit code
 */
int main(int argc, char* argv[]){
  // Checking input format
  if(argc < 2){
    printf("Please enter rna sequence after the file name.\n");
    return 0;
  }
  // For printing of ouput in correct order
  LARGE_INTEGER frequency; /**< Store ticks per second*/
  LARGE_INTEGER t1, t2;  /**< Stores the number of ticks*/
  double elapsed_time; /**< Stores the net elapsed time*/

  // get ticks per second
  QueryPerformanceFrequency(&frequency);
  // start timer
  QueryPerformanceCounter(&t1);

  setbuf(stdout, NULL);
  char *rna = argv[1]; /**< The input RNA sequence.*/
  int n; /**< Stores length of the input RNA sequence.*/
  int m; /**< Stores the number of pairs in the output.*/

  // Preprocessing input
  lrtrim(rna);
  n = strlen(rna);

  char* out = (char*)calloc(n+1, sizeof(char)); /**< Stores the output for Varna visualization*/

  // Initializing out array
  for(int i = 0; i < n; i++){
    out[i] = '.';
  }
  out[n] = 0;

  // Checking base condition
  if(n > 5){
    int** to_clear = (int**)calloc(n*n, sizeof(int*)); /**< Stores pointers of base pairs to be freed*/
    int to_clear_i = 0; /**< Index of the next entry in to_clear array*/
    int**** OPTa = (int****)calloc(n-5, sizeof(int***)); /**< Dynamic Programming Array, stores the pairs maximizing the maximum number of pairs found till now*/
    int** OPT = (int**)calloc(n-5, sizeof(int*)); /**<Dynamic Programming Array, stores maximum number of pairs found till now*/
    
    // Initializing OPT and OPTa
    for(int i = 0; i < n-5; i++){
      OPT[i] = (int*)calloc(n-5, sizeof(int));
      OPTa[i] = (int***)calloc(n-5, sizeof(int**));
      for(int j = 0; j < n-5; j++){
        OPTa[i][j] = (int**)calloc(n, sizeof(int*));
      }
    }

    // run the algorithm on the input
    rnaFold(to_clear, to_clear_i, OPTa, OPT, n, rna);

    // Printing the answers
    printf("Answer: %d Pairs Matched\n", OPT[0][n-6]);
    m = OPT[0][n-6];
    printf("::Pairs matched::\n");

    // Building the output array
    buildOut(m, OPTa, out, n);

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
  char command[1000]; /**< Stores the command to run the Varna tool*/
  sprintf(command, "java -cp assets/bin/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN %s -o viz.svg -algorithm naview", rna, out);
  int result = system(command); /**< Stores the result of running the command*/
  printf("Vizualization saved in the output file\n");

  // free output array
  free(out);

  // stop timer
  QueryPerformanceCounter(&t2);
  // compute and print the elapsed time in millisec
  elapsed_time = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
  printf("Number of Bases: %d, Number of pairs made: %d\n", n, m);
  printf("Output computed in: %f ms.\n", elapsed_time);
  return 0;
}