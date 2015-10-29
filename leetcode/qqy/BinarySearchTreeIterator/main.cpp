/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class BSTIterator {
    stack<TreeNode*> s;
    TreeNode* node = NULL;
public:
    BSTIterator(TreeNode *root) {
        while(root){
            s.push(root);
            root = root->left;
        }
        node = root;
    }

    /** @return whether we have a next smallest number */
    bool hasNext() {
         return (!s.empty() || (node != NULL));
    }

    /** @return the next smallest number */
    int next() {
        while(node){
            s.push(node);
            node = node->left;
        }
        node = s.top();
        s.pop();
        int res = node->val;
        node = node->right;
        return res;
    }
};

/**
 * Your BSTIterator will be called like this:
 * BSTIterator i = BSTIterator(root);
 * while (i.hasNext()) cout << i.next();
 */

//The iteratoer is an object that captures the traversal info
