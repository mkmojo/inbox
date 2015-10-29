/**
 * Definition for binary tree with next pointer.
 * struct TreeLinkNode {
 *  int val;
 *  TreeLinkNode *left, *right, *next;
 *  TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
 * };
 */

//This solution can handle non-perfect binary tree as well
class Solution {
public:
    void connect(TreeLinkNode *root) {
        if(root == NULL) return ;
        queue<TreeLinkNode*> q;
        q.push(root);
        q.push(NULL);
        TreeLinkNode * lastNode = NULL;
        while(!q.empty()){
            TreeLinkNode *cnt_node = q.front();
            q.pop();
            //pay attention to termination condition
            //only when q is not empty, should you push in specail node, 
            //otherwise will not terminate
            if(cnt_node == NULL && !q.empty()) q.push(cnt_node);
            if(cnt_node && cnt_node->left)  q.push(cnt_node->left);
            if(cnt_node && cnt_node->right) q.push(cnt_node->right);
            
            if(lastNode) lastNode->next = cnt_node;
            lastNode = cnt_node;
        }
    }
};
