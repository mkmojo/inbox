#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

class Solution {
    public:
        string addBinary(string a, string b) {
            string l = a.size() > b.size() ? a : b;
            string s = a.size() <= b.size() ? a : b;
            reverse(l.begin(),l.end());
            reverse(s.begin(),s.end());

            string res="";
            _carry = '0';
            for(int i=0; i<l.size();i++){
                if(i<s.size()){
                    res.append(string(1, add(s[i],l[i])));
                    continue;
                }
                res.append(string(1, add('0', l[i])));
            }
            res += (_carry == '0' ? "": "1");
            reverse(res.begin(), res.end());
            return res;
        }

    private:
        char _carry;
        char add(char a, char b){
            int x = a - '0';
            int y = b - '0';
            int c = _carry - '0';

            char res = (x + y + c) % 2 + '0';
            _carry = (x + y + c) / 2 ? '1' : '0';
            return res;
        }
};

int main()
{
    Solution s;
    cout << "test 10, 1000:  " << (s.addBinary("10", "1000") == "1010"? "pass!" : "fail!") << endl;
    cout << "test 1, 1:  " << (s.addBinary("1", "1") == "10"? "pass!" : "fail!") << endl;
    return 0;
}
