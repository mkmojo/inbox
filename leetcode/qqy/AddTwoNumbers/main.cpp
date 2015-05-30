#include <iostream>
#include <vector>
using namespace std;

struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};

ListNode* makeList(vector<int> vals){
    ListNode* head = NULL;
    ListNode* p = NULL;

    for( int i=0;i<vals.size();i++){
        ListNode* node = new ListNode(vals[i]);
        if(i==0){ 
            head = node;
            p = head;
        }
        p->next = node;
        p = p->next;
    }
    return head;
}

class Solution {
    //calculate result and set value to r
    void add(ListNode *x, ListNode *y, int *carry, ListNode *r){
        int res = (x?x->val:0) + (y?y->val:0) + *carry;
        if(res >= 10){
            r->val = res%10;
            *carry = 1;
        }else{
            r->val = res;
            *carry = 0;
        }
        cout << "add-> res: " << res << " r-val: " << r->val << " carry: " << *carry<<endl;
    }
    public:
    void printList(ListNode *head, string listName){
        ListNode *p = head;
        cout << listName << ": ";
        while(p != NULL){
            cout << p->val <<" ";
            p = p->next;
        }
        cout << endl;
    }
    ListNode* addTwoNumbers(ListNode* l1, ListNode* l2) {
        if(l1 == NULL) return l2;
        if(l2 == NULL) return l1;
        ListNode *l3 = new ListNode(0);
        ListNode *r = l3;
        int carry = 0;

        while(l1 && l2){
            add(l1, l2, &carry, r);
            printList(l1, "l1");
            printList(l2, "l2");
            printList(l3, "l3");
            r->next = new ListNode(-1);
            r = r->next;
            l1 = l1->next;
            l2 = l2->next;
        }

        while(carry == 1){
            if(!l1 && !l2 ){
                r->next = new ListNode(carry);
                carry = 0;
            }else if(!l1 && l2){
                add(l1, l2, &carry, r);
                l2=l2->next;
                r->next = new ListNode(-1);
                r = r->next;
            }else if(l1 && !l2){
                add(l1, l2, &carry, r);
                l1=l1->next;
                r->next = new ListNode(-1);
                r = r->next;
            }else{
                cout << "error" <<endl;
            }
        }
        return l3;
    }
};

int main()
{
    Solution s;
    vector<int> x = {9, 9, 9};
    vector<int> y = {1, 0, 0};
    ListNode* p = makeList(x);
    ListNode* q = makeList(y);
    s.printList(s.addTwoNumbers(p,q), "ans");

    return 0;
}
