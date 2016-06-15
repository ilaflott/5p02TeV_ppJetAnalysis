#include <cstdlib>
#include <ctime>
#include <iostream>


using namespace std;

int main(int argc, char *argv[]){
  //cout<<"hello world, testing vector stuff. don't break!" << endl;
  int timeint=time(NULL);
  cout << "timeint="<<timeint<<endl;
  srand(timeint);

  for (int i=0;i<10;i++)
    {
      cout<<"rand="<<rand()<<endl;
    }  
  return 0;
}
