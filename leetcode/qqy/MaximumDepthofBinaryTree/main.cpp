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

    void do_maxDepth(TreeNode* root, int level, int& cnt_max){
        if(level > cnt_max) cnt_max = level;
        if(root->left != NULL) do_maxDepth(root->left, level + 1, cnt_max);
        if(root->right != NULL) do_maxDepth(root->right, level + 1, cnt_max);
    }
public:
    int maxDepth(TreeNode* root) {
        if (root == NULL) return 0;
        int maxLevel = 0;
        do_maxDepth(root, 1, maxLevel);
        return maxLevel;
    }
};
