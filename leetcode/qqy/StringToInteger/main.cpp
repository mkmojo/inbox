#include<string>
#include<iostream>
#include<stdio.h>
using namespace std;

//#define MAX_INT  2147483647
//#define MIN_INT  (-2147483647 - 1)
class Solution {
    bool isSpace(char x){
        static const string whiteSpaces=" \t\n\r";
        for(int i=0;i<whiteSpaces.size();i++){
            if(x == whiteSpaces[i]) return true;
        }
        return false;
    }
    bool isNumber(char x){
        return x>='0' && x<='9';
    }
    public:
    const int MAX_INT=(1<<(sizeof(int)*8 - 1)) - 1;
    const int MIN_INT=(1<<(sizeof(int)*8 - 1));
    int myAtoi(string str) {
        int i=0;
        while(isSpace(str[i++]));
        str = str.substr(i-1, str.size() - i + 1);
        //cout << "DEBUG: str: " << str<<endl;
        bool isNeg = false;
        if(str[0] == '-' || str[0] == '+'){
            isNeg = (str[0] == '-');
            //only when the first charactor is a sign
            str = str.substr(1,str.size()-1);
        }

        int sum = 0;
        for(int i=0;isNumber(str[i]);i++){
            int d = str[i] - '0';
            if(isNeg){
                if(-sum < (MIN_INT + d)/10) {
                    //cout << "return MIN_INT ";
                    return MIN_INT;
                }
            }else{
                if(sum > (MAX_INT - d)/10) {
                    //cout << "return MAX_INT ";
                    return MAX_INT;
                }
            }
            sum  = sum*10 + d;
            //cout << "DEBUG: " << sum <<endl;
        }
        return isNeg ? -sum : sum; 
    }
};

int main()
{
    Solution s;
    cout << "MAX_INT " << s.MAX_INT <<endl;
    cout << "MIN_INT " << s.MIN_INT <<endl;
    printf("%s = %d\n", "'123'",s.myAtoi("123"));
    printf("%s = %d\n", "'2147483648'",s.myAtoi("2147483648"));
    printf("%s = %d\n", "'-2147483648'",s.myAtoi("-2147483648"));
    //cout << "'123' ="  << s.myAtoi("123") <<endl;
    //cout << "'2147483648' "<< "=" << s.myAtoi("2147483648")<<endl;
    //cout << "'-2147483648' " << "="<< s.myAtoi("-2147483648")<<endl;
    return 0;
}
