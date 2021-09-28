//
// Created by Jan Schäppi on 22.09.21.
//

//#include <vector>
using namespace std;

template <typename T>
int maxOfDifference(const vector<T>& v1, const vector<T>& v2){
    if(v1.size() != v2.size()) cout << "Error: v1 must have the same size as v2" << endl;
    if(v1.size() == 0) cout << "Length of v1 is zero! This is not allowed" << endl;
    int maximum = v1[0]-v2[0];
    for(int i = 1; i < v1.size(); ++i){
        if(v1[i]-v2[i] > maximum) maximum = v1[i]-v2[i];
    }
    return maximum;
}

template <typename T>
vector<T> strictLowerDiagonalMatrix(const int n, T init){
    vector<T> vec((n*(n-1))/2, init);
    return vec;
}

int positionInStrictLowerDiagonalMatrix(const int i, const int j){
    return (i*(i-1))/2 + j;
}

template <typename T>
void printStrictLowerDiagonalMatrix(const vector<T>& matrix_mxn, const int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(j < i) cout << matrix_mxn[positionInStrictLowerDiagonalMatrix(i,j)] << " ";
            else cout << "?" << " ";
        }
        cout << endl;
    }
}


bool endswith(string str, string substring){
    if(str.size() < substring.size()) return false;
    return str.substr(str.size()-substring.size(),substring.size()) == substring;
}

bool equals(const char* s1, const char* s2, const int l){
    for(int i = 0; i < l; ++i){
        if(*(s1++) != *(s2++)) return false;
    }
    return true;
}

string to_string(string str){
    return str;
}

template<typename T>
void print_vectors(vector<vector<T>> vecs){

    // Figure out the maximum length
    int size = vecs[0].size();
    for(int i = 0; i < vecs.size(); ++i){
        if(vecs[i].size() != size){
            cout << "ERROR while print_vectors: Vectors have to be of the same size." << endl;
            return;
        }
    }

    string printout; int max_chars = int(100/vecs.size())-3;
    for(int i = 0; i < size; ++i){
        cout << i << " : ";
        for(int j = 0; j < vecs.size(); ++j){
            printout = to_string(vecs[j][i]);
            if(printout.size() > 100) cout << printout.substr(0,max_chars) << "...";
            else cout << printout;
            if(j < vecs.size()-1) cout << " | ";
        }
        cout << endl;
    }
}

template<typename T>
void print_hotizontal(vector<T>& vec){
    for(int j = 0; j < vec.size(); ++j) cout << vec[j] << " ";
    cout << endl;
}


template<typename T>
void print_vectors(vector<T>& vec1, vector<T>& vec2, vector<T>& vec3){
    vector<vector<T>> vecs;
    vecs.push_back(vec1); vecs.push_back(vec2); vecs.push_back(vec3);
    return print_vectors(vecs);
}

template<typename T>
void print_vectors(vector<T>& vec1, vector<T>& vec2){
    vector<vector<T>> vecs;
    vecs.push_back(vec1); vecs.push_back(vec2);
    return print_vectors(vecs);
}


template<typename T>
void print_vectors(vector<T>& vec1){
    vector<vector<T>> vecs;
    vecs.push_back(vec1);
    return print_vectors(vecs);
}