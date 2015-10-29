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
    vector<int> linearize(TreeNode* node){
        stack<TreeNode*> s;
        vector<int> ans;
        while(!s.empty() || node != NULL){
            if(node != NULL){
                s.push(node);
                node = node->left;
            }else{
                node = s.top();
                s.pop();
                ans.push_back(node->val);
                node = node->right;
            }
        }
        return ans;
    }
public:
    int kthSmallest(TreeNode* root, int k) {
        vector<int> arr = linearize(root);
        return arr[k-1];
    }
};

//My idea is to traverse the tree first to get an iteration
//Then directly index into it to get the element. 
//Do not yet know how to do it recursivlly. 
//The challenge for me is that I do not know 

// Here is the recursive version
// The idea is to throw out half of elements everytime recursively calls
// itself.
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
    int treeSize(TreeNode* node){
        if(node == NULL) return 0;
        int nl = treeSize(node->left);
        int nr = treeSize(node->right);
        return nl + nr + 1;
    }
public:
    int kthSmallest(TreeNode* root, int k) {
        int nl = treeSize(root->left);
        if(k == nl + 1) return root->val;
        else if(k > nl + 1) return kthSmallest(root->right, k-nl-1);
        else return kthSmallest(root->left, k);
    }
};

//ref answer is neat in that  
//whenever the program visit some node, that means I have resolve
//a small part of my problem and hense the size of k reduces by 1.
//This way the program does not need to save all results. 

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
    int kthSmallest(TreeNode* root, int k) {
        stack<TreeNode*> ps;
        while(!ps.empty() || root != NULL){
            while(root){
                ps.push(root);
                root = root->left;
            }
            root = ps.top();
            ps.pop();
            k--;
            if(k == 0) return root->val;
            root = root->right;
        }
    }
};
