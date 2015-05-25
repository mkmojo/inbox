#include<iostream>
#include<assert.h>
#include<string>
#include<algorithm>
#include<vector>
using namespace std;
class Solution {
    void do_reverse(string &s){
        int i = s[0]=='-' ? 1 : 0;
        int j = s.size() - 1;
        while(i<j){
            char c = s[i];
            s[i] = s[j];
            s[j] = c;
            i++;
            j--;
        }
    }
    //takes in a string that represnts positive number only
    int str_2_long(string s){
        int sum = 0;
        const int MAXINT = (1L<<(sizeof(int)*8-1))-1;
        bool isNeg = s[0] == '-';
        for(int i = isNeg ? 1 : 0; i<s.size();i++){
            char c = s[i];
            if(sum> MAXINT/10) return 0; // <--- I added this line
            sum = (c - '0') + sum*10;
        }
        
        return isNeg ? -sum : sum;
    }
public:
    int reverse(int x) {
        string s = to_string(x);
        do_reverse(s);
        return str_2_long(s);
    }
};


int main()
{
    Solution s;
    int i = 1534236469; 
    assert(s.reverse(123) == 321);
    assert(s.reverse(-321) == -123);
    assert(s.reverse(1534236469) == 0 );
    cout<< "test all passed" <<endl;
    return 0; 
}
