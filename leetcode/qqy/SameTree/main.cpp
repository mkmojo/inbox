/**
 * Definition for a binary tree node.
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    bool isSameTree(TreeNode* p, TreeNode* q) {
        if(p == NULL && q == NULL) return true;
        // Remember to check for p and q before -> 
        // into p and q node. Otherwise Runtime Error will occur
        if(p && q && p->val == q->val){
            return isSameTree(p->left, q->left) &&
            isSameTree(p->right, q->right);
        }
        return false;
    }
};
