nclude <iostream>
#include <string>

using namespace std;

class Solution {
    int reverse(int x){
        int res = 0;
        while(x != 0){
            int r = x%10;
            res = res*10 + r;
            x /= 10;
        }
        return res;
    }
public:
    bool isPalindrome(int x) {
        return x>=0 && reverse(x) == x;
    }
};

int main()
{
    Solution s;
    cout << (s.isPalindrome(-2147447412) == false ? "pass" : "fail") <<endl;
}

