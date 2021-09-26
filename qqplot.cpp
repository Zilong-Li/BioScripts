#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

string timestamp()
{
    time_t t = time(NULL);
    char *s = asctime(localtime(&t));
    s[strlen(s) - 1] = '\0';
    string str(s);
    str = string("[") + str + string("] ");
    return str;
}

int main(int argc, char *argv[])
{
    cerr << timestamp() << "start processing.\n";
    int bins = 1000;
    double cutoff = 1e-4;
    double c = -log10(cutoff);
    int64_t N = 0, i = 0;
    auto allbins = new double[bins][2];
    for (i = 0; i < bins; i++) {allbins[i][0] = -1;allbins[i][1] = 0;}
    double p;
    vector<double> sigp;
    ios_base::sync_with_stdio(false);
    while ( cin >> p) {
        N += 1;
        if (N % 100000000 == 0) { cerr << timestamp() << "reached the " << N << " line.\n"; }
        p = -log10(p);
        if (p > c) {
            sigp.push_back(p);
        } else {
            i = bins - int(floor((double) bins * p / c)) - 1;
            i == -1 ? i = 0 : i;
            allbins[i][0] = fmax(allbins[i][0], p);
            allbins[i][1] += 1;
        }
    }
    cerr << timestamp() << "the total number of lines is " << N << ".\n";

    int64_t size = sigp.size();
    sort( sigp.begin(), sigp.end(), std::greater<double>() );
    vector<double> expp;
    for (i = 0; i < size; ++i) {
        expp.push_back(-log10( (double)(i + 1 - 0.5) / N ));
    }

    for (i = 0; i < bins; ++i){
        if (allbins[i][0] != -1){
            sigp.push_back(allbins[i][0]);
            expp.push_back(-log10( (double)(allbins[i][1] + size - 0.5) / N) );
            size += allbins[i][1];
        }
    }

    cout.precision(17);
    for (i = 0; i < sigp.size(); ++i){
        cout << sigp[i] << "," << expp[i] << "\n";
    }
    
    cerr << timestamp() << "program finished.\n";
    delete[] allbins;
    
    return 0;
}