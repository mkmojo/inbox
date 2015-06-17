#include <iostream>
using namespace std;
typedef struct TreeNode{
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode(int x) : val(x), left(NULL), right(NULL) {}
}TreeNode;

class Solution {
    void do_invert(TreeNode *root, TreeNode *newRoot){
        if(root == NULL) return;
        if(root->left){
            TreeNode* newNode = new TreeNode(root->left->val);
            newRoot->right = newNode;
            do_invert(root->left, newRoot->right);
        }
        if(root->right){
            TreeNode* newNode = new TreeNode(root->right->val);
            newRoot->left = newNode;
            do_invert(root->right, newRoot->left);
        }
    }
    public:
    TreeNode* invertTree(TreeNode* root) {
        TreeNode *newRoot = new TreeNode(root->val);
        do_invert(root, newRoot);
        return newRoot;
    }
};

int main()
{
    Solution s;
    return 0;
}
