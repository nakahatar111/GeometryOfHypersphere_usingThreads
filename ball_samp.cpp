#include <iostream>
#include <stdio.h>
#include <chrono>
#include <vector>
#include <math.h>
#include <omp.h> 
#include <random>
#include <iomanip>
using namespace std;
using namespace std::chrono;

vector<vector<int>> ball_samp_new_type(int dim, int numPoints){
  default_random_engine generator;
  vector<vector<int>> resultPerDim(dim-1, vector<int>(100, 0)); 
  for(int d = 2; d <= dim; d++){
    for(int p = 0; p < numPoints; p++){
      double dist = 0.0;

      do{
        double coordSquaredSum = 0.0;
        double range = 1.0;
        for(int numCoord = 0; numCoord < d; numCoord++){
          range = range - coordSquaredSum;
          cout << "range: +-" << range << endl;
          uniform_real_distribution<double> distribution(-1*sqrt(range),sqrt(range));  
          coordSquaredSum += pow(distribution(generator),2);
        }
        dist = sqrt(coordSquaredSum);
        cout << "here: " << coordSquaredSum << endl;
      }while (dist >= 1);
      cout << "dist:  " << dist << endl; 
      resultPerDim.at(d-2).at((int)(dist*100)) += 1;
    }
  }
  return resultPerDim;
}

vector<vector<int>> ball_samp_new(int dim, int numPoints){
  default_random_engine generator;
  uniform_real_distribution<double> distribution(-1.0,1.0);
  vector<vector<int>> resultPerDim(dim-1, vector<int>(100, 0)); 
  for(int d = 2; d <= dim; d++){
    for(int p = 0; p < numPoints; p++){
      int numCoord = 0;
      double coordSquaredSum = 0.0;

      while (numCoord < d){
        double randNum = pow(distribution(generator),2);
        if(sqrt(coordSquaredSum + randNum) < 1){
          coordSquaredSum = coordSquaredSum + randNum;
          numCoord++;
        }
      }

      resultPerDim.at(d-2).at((int)(sqrt(coordSquaredSum)*100)) += 1;
    }
  }
  return resultPerDim;
}

vector<vector<int>> ball_samp_new_parallel(int dim, int numPoints){
  default_random_engine generator;
  uniform_real_distribution<double> distribution(-1.0,1.0);
  vector<vector<int>> resultPerDim(dim-1, vector<int>(100, 0));
  #pragma omp parallel private(generator, distribution)
  {
    generator.seed(omp_get_thread_num());
    #pragma omp for
    for(int d = 2; d <= dim; d++){
      for(int p = 0; p < numPoints; p++){
        int numCoord = 0;
        double coordSquaredSum = 0.0;

        while (numCoord < d){
          double randNum = pow(distribution(generator),2);
          if(sqrt(coordSquaredSum + randNum) < 1){
            coordSquaredSum = coordSquaredSum + randNum;
            numCoord++;
          }
        }

        resultPerDim.at(d-2).at((int)(sqrt(coordSquaredSum)*100)) += 1;
      }
    }
  }
  return resultPerDim;
}


vector<vector<int>> ball_samp_seq(int dim, int numPoints){
  default_random_engine generator;
  uniform_real_distribution<double> distribution(-1.0,1.0);
  vector<vector<int>> resultPerDim(dim-1, vector<int>(100, 0)); 
  for(int d = 2; d <= dim; d++){
    for(int p = 0; p < numPoints; p++){
      double dist;

      do{
        dist = 0.0;
        for(int numCoord = 0; numCoord < d; numCoord++)
          dist += pow(distribution(generator),2);
        dist = sqrt(dist);
      }while (dist >= 1);

      
      resultPerDim.at(d-2).at((int)(dist*100)) += 1;
    }
  }
  return resultPerDim;
}

vector<vector<int>> ball_samp_parallel(int dim, int numPoints){
  default_random_engine generator;
  uniform_real_distribution<double> distribution(-1.0,1.0);
  vector<vector<int>> resultPerDim(dim-1, vector<int>(100, 0)); 
  #pragma omp parallel private(generator)
  {
    generator.seed(omp_get_thread_num());
    #pragma omp for collapse(2) schedule(dynamic, 4)
    for(int d = 2; d <= dim; d++){ 
      for(int p = 0; p < numPoints; p++){
        double dist;

        do{
          dist = 0.0;
          for(int numCoord = 0; numCoord < d; numCoord++)
            dist += pow(distribution(generator),2);
          dist = sqrt(dist);
        }while (dist >= 1);
        
        #pragma omp atomic
        resultPerDim.at(d-2).at((int)(dist*100)) += 1;
      }
    }
  }
  return resultPerDim;
}

void test_fn_seq(){
  int sum = 0;
  for(int d = 2; d<=16;d++){
    for(int p = 0; p < 100; p++){
      int i = 0;
      for(int numCoord = 0; numCoord < 10000; numCoord++){
        i= numCoord+i;
      }
      sum += i^2;
    }
  }
  cout << sum << endl;
}
void test_fn_parallel(){
  int sum = 0;
  #pragma omp parallel reduction(+:sum)
  {
    cout << sum;
    #pragma omp for
    for(int d = 2; d<=16;d++){
      for(int p = 0; p < 100; p++){
        int i = 0;
        for(int numCoord = 0; numCoord < 10000; numCoord++){
          i= numCoord+i;
        }
        sum += i^2;
      }
    }
  }
  cout << sum << endl;
}


int main (int argc, char *argv[]){

  int dimension = 16;
  int numpoints = 100;
  /*
  auto start7 = high_resolution_clock::now();
  //ball_samp_new_type(dimension, numpoints);
  auto stop7 = high_resolution_clock::now();
  auto duration7 = duration_cast<microseconds>(stop7 - start7);
  cout << "Test Time: " << (double)duration7.count() << endl;

  auto start4 = high_resolution_clock::now();
  test_fn_seq();
  auto stop4 = high_resolution_clock::now();
  auto duration4 = duration_cast<microseconds>(stop4 - start4);
  cout << "Test1 Time: " << (double)duration4.count() << endl;

  auto start5 = high_resolution_clock::now();
  test_fn_parallel();
  auto stop5 = high_resolution_clock::now();
  auto duration5 = duration_cast<microseconds>(stop5 - start5);
  cout << "Test2 Time: " << (double)duration5.count() << endl;
  */
  
  auto start1 = high_resolution_clock::now();
  vector<vector<int>> output1 = ball_samp_seq(dimension, numpoints); 
  auto stop1 = high_resolution_clock::now();
  auto duration1 = duration_cast<microseconds>(stop1 - start1);
  
  auto start2 = high_resolution_clock::now();
  vector<vector<int>> output2 = ball_samp_parallel(dimension, numpoints);
  auto stop2 = high_resolution_clock::now();
  auto duration2 = duration_cast<microseconds>(stop2 - start2);

  cout << fixed << setprecision(6);

  cout << "Output 1 Time: " << (double)duration1.count()/1000000 << "s" << endl;
  cout << "Output 2 Time: " << (double)duration2.count()/1000000 << "s" << endl;
  cout << fixed << setprecision(2);
  cout << "SpeedUp: x" << (double)duration1.count()/duration2.count() << endl << endl;
  
  cout << "OUTPUT 1 (Sequential)"<< endl;
  for(int j = 0; j < output1.size(); j++){
    cout << "Dimension: " << j+2 << endl;
    int sum = 0;
    for(int i = 0; i < 100; i++){
      cout << "  At " << (double)(i/100.0) << " count:";
      int count = output1.at(j).at(i);
      for(int a = 0; a < count; a++)
        cout << "-";
      if (count == 0)
        cout << " " << count << endl;
      else
        cout << " (" << count << "/" << numpoints << ")" << endl;
      sum += output1.at(j).at(i);
    }
    cout << "Num Points: " << sum << endl;
  }
  
  cout << endl << endl << endl << "OUTPUT 2 (Parallel)"<< endl;
  for(int j = 0; j < output2.size(); j++){
    cout << "Dimension: " << j+2 << endl;
    int sum = 0;
    for(int i = 0; i < 100; i++){
      cout << "  At " << (double)(i/100.0) << " count:";
      int count = output2.at(j).at(i);
      for(int a = 0; a < count; a++)
        cout << "-";
      if (count == 0)
        cout << " " << count << endl;
      else
        cout << " (" << count << "/" << numpoints << ")" << endl;
      sum += output2.at(j).at(i);
    }
    cout << "Num Points: " << sum << endl;
  }
}