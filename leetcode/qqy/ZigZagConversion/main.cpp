#include<iostream>
#include<string>
#include<vector>
#include<stdio.h>
using namespace std;
class Solution {
    public:
        string convert(string s, int numRows) {
            if(s.size() == 0) return "";
            if(numRows == 1) return s;

            int x = 0;
            int down = 0;
            vector<string> v(numRows);
            for(int i=0;i<s.size();i++){
               if(i%((numRows-1)*2)< numRows - 1 ){ 
                   //down
                    v[x]+=s[i];
                    x+=1;
                } else{ 
                    //up-right
                    v[x]+=s[i];
                    x-=1;
                }
            }
            string ans = "";
            for(int i=0; i<numRows;i++){
                ans += v[i];
            }
            return ans;
        }
};
int main()
{
    Solution s;
    cout << s.convert("AB", 1) <<endl;
    cout << s.convert("ABCDEFG", 2) <<endl;
    cout << s.convert("ABCDEFG", 3) <<endl;

    return 0;
}
