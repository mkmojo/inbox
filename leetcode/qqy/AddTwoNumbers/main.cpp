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
            continue;
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
        cout << "DEBUG: add-> res: " << res << " r-val: " << r->val 
            << " carry: " << *carry<<endl;
    }

    int readAndAdvance(ListNode * &l){
        int res = 0;
        if(l){
            res = l->val;
            l = l->next;
        }
        return res;
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
            //build node
            r->next = new ListNode(-1);
            r = r->next;
            add(l1, l2, &carry, r);
            //printList(l1, "l1");
            //printList(l2, "l2");
            //printList(l3->next, "l3");
            l1 = l1->next;
            l2 = l2->next;
        }

        while(carry == 1){
            r->next = new ListNode(-1);
            r = r->next;
            if(!l1 && !l2 ){
                r->val = carry;
                carry = 0;
            }else if(!l1 && l2){
                add(l1, l2, &carry, r);
                l2=l2->next;
            }else if(l1 && !l2){
                add(l1, l2, &carry, r);
                l1=l1->next;
            }else{
                cout << "error" <<endl;
            }
        }

        l1 ? r->next = l1 : l2 ? r->next = l2 : NULL;
        return l3->next;
    }

    ListNode* addTwoNumbers1(ListNode* l1, ListNode* l2){
        ListNode *h = NULL, **pp = &h;
        int carry = 0, x = 0, y = 0, sum = 0;

        while(l1 || l2){
            x = readAndAdvance(l1);
            y = readAndAdvance(l2);

            sum = x + y + carry;

            ListNode* node = new ListNode(sum%10);
            *pp = node; 
            pp = &(node->next);

            carry = sum / 10;
        }

        if(carry == 1){
            ListNode* node = new ListNode(1);
            *pp = node; 
            pp = &(node->next);
        }
        return h;
    }
};

int main()
{
    Solution s;
    vector<int> x = {2, 4, 3, 3,4,5,6,7,2,5,2,2,3};
    vector<int> y = {5, 6, 4};
    vector<int> m = {9, 9, 9};
    vector<int> n = {1};
    vector<int> a = {0};
    vector<int> b = {0};
    ListNode* p = makeList(m);
    ListNode* q = makeList(n);
    s.printList(s.addTwoNumbers1(p,q), "ans");

    return 0;
}
