#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]){
  //cout<<"hello world, testing vector stuff. don't break!" << endl;
  vector<int> intVector1 (10,0);//10 ints value 0
  int myints[]={3,4,5,6,7};
  vector<int> intVector2 (myints, myints+sizeof(myints)/sizeof(int) );//declare w/ array
  
  vector<int> *intVectorPoint;

  //cout<<"intVector2="<<intVector2 <<endl;//doesn't work! vectors are templates/structures so they have no sensible cout 
  cout<<  "&intVector2    ="<<&(intVector2) <<endl;//does work!
  intVectorPoint=&(intVector2);
  cout << "intVectorPoint ="<<intVectorPoint<<endl;//i am address of intVector2
  cout << "&intVectorPoint="<<&intVectorPoint<<endl;//i am address of intVectorPoint
  //cout << "*intVectorPoint="<<*intVectorPoint<<endl;//doesn't work! vectors are templates/structures so they have no sensible cout 

  cout << " myints="<<myints<<endl; //i am the address of the start of the array
  cout << "&myints="<<&myints<<endl;//so am i
  cout << "*myints="<<*myints<<endl;//i am the contents of the first element in the array
  cout << "*myints="<<*(myints+1)<<endl;//i am the contents of the second element in the array
  
  
  return 0;
}
