#include <stdio.h>
#include <stdlib.h>

int main(){
  FILE *in = fopen("file_names.dat", "r");
  FILE *out = fopen("run_files.sh","w");
  char  c;
  int n=0;
  while((c=getc(in))!=EOF){
    ungetc(c,in);
    n=0;
    char filename[100]="";
    while ((c=getc(in))!='\n' && c!= EOF){
      filename[n]=c;
      n++;
    }
    fprintf(out, "%s%s \n", "./a.out ", filename);
  }
  fclose(in);
  fclose(out);
  return 0;
}
