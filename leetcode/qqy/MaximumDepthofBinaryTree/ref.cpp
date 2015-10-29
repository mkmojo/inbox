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
    int maxDepth(TreeNode* root) {
        if (root == NULL) return 0;
        if (!root->left && !root->right){
            return 1;
        }
        // Wrong int d_left = maxDepth(root->left); 
        // Wrong int d_right = maxDepth(root->right);
        int d_left = 1, d_right = 1;
        d_left  = maxDepth(root->left);
        d_right  = maxDepth(root->right);

        return d_left > r_left ? d_left : d_right;
    }
};
