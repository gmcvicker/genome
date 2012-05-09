
import sys


class TreeNode(object):
    def __init__(self):
        self.name = None
        self.next = None
        self.branch = None
        self.branch_len = None
        
    def is_leaf(self):
        """Returns True if this Node is a leaf"""
        return (self.next == self)

    def is_root(self):
        """Returns True if this is the root node"""
        return (self.branch is None)


    def find_node(self, name):
        """Recursively searches descendents of this tree for a node with the
        given name. Returns None if no node matches the name exactly. If
        multiple nodes match, only the first matching node is returned"""
        if (self.name is not None) and self.name == name:
            return self

        cur = self.next
        while cur != self:
            if (cur.name is not None) and (cur.name == name):
                return cur
            else:
                # recursively search descendents
                node = cur.branch.find_node(name)

                if node is not None:
                    return node

            cur = cur.next

        return None


    def find_root(self, from_branch=None):
        """Returns the root of the provided tree or None if no root is
        found"""
        if self.is_root():
            sys.stderr.write("root")
            return self

        if not from_branch:
            res = self.branch.find_root(from_branch=self)
            if res:
                return res

        cur = self.next
        while cur != self:
            if cur.is_root():
                return cur

            res = cur.branch.find_root(from_branch=cur)
            if res:
                return res
            cur = cur.next

        return None
            


    def copy(self, from_branch=None):
        """Make a copy of this tree"""

        # copy this node
        new_node = TreeNode()
        new_node.name = self.name

        if self.branch:
            new_node.branch_len = self.branch_len
            if from_branch:
                new_node.branch = from_branch
            else:
                new_node.branch = self.branch.copy(from_branch=new_node)

        # copy linked nodes
        cur = self.next
        new_prev = new_node
        new_cur = None
        
        while cur != self:
            new_cur = TreeNode()
            new_cur.name = cur.name
            if cur.branch:
                # recursively copy branches off of each node
                new_cur.branch_len = cur.branch_len
                new_cur.branch = cur.branch.copy(from_branch=new_cur)
            new_prev.next = new_cur
            new_prev = new_cur
            cur = cur.next

        # join up circle of linked nodes
        new_prev.next = new_node
        return new_node
            


    def remove(self):
        """Returns a copy of this tree, with this node removed. The returned
        reference is the root node of the new tree, or an arbitrary node
        if there is no root."""
        node = self.copy()
        
        if node.is_leaf():
            # remove parent
            if node.branch:
                node = node.branch
            else:
                # this tree was only a leaf, and it is now removed
                return None
        
        prev = None
        cur = node.next
        n = 1
        # find previous node
        while cur != node:
            if cur.is_leaf():
                raise ValueError("internal leaf?")
            prev = cur
            cur = cur.next
            n += 1
        
        # remove this node from next list
        prev.next = node.next

        if n == 3:
            # there are only two nodes left, we need to
            # merge these to make a longer branch
            parent = node.next.branch
            child = prev.branch

            if child is None and parent is None:
                raise ValueError("Tree defines multiple roots")

            if child is None:
                # parent is now root of tree
                parent.branch = None
                parent.branch_len = None
                return parent
            elif parent is None:
                # child is now root of tree
                child.branch = None
                child.branch_len = None
                return child

            # join branches
            parent.branch = child
            child.branch = parent

            if child.branch_len and parent.branch_len:
                child.branch_len = parent.branch_len + child.branch_len
                parent.branch_len = child.branch_len
            else:
                child.branch_len = None
                parent.branch_len = None

            root = parent.find_root()
            if root:
                return root
            else:
                return parent

        root = prev.find_root()
        if root:
            return root
        else:
            return prev

        
    
    def remove_leaf(self, name):
        """Returns a copy of this tree, with the leaf of the provided name
        removed. The returned node is the root of the tree, unless the tree
        is unrooted, in which case an arbitrary node is returned.
        Raises a ValueError if a node with the provided name can not be
        found, or if the named node is not a leaf."""
        leaf = self.find_node(name)

        if leaf is None:
            raise ValueError("no node of name %s exists" % name)

        if not leaf.is_leaf():
            raise ValueError("node %s is not a leaf" % name)

        return leaf.remove()


            
    def remove_leaves(self, names):
        """Returns a copy of this tree with the leaves matching names in
        the provided list removed"""
        if len(names) == 0:
            return self.copy()

        tree = self
        for name in names:
            tree = tree.remove_leaf(name)

        return tree
        

        

    def __str__(self):
        """Converts tree to a newick string representation"""
        tree_str = ""
        if self.branch_len is None:
            branch_len_str = ""
        else:
            branch_len_str = ":%g" % self.branch_len

        if self.is_leaf():
            # this is a leaf
            if self.name is None:
                name_str = ""
            else:
                name_str = self.name
            tree_str = name_str + branch_len_str
        else:
            # internal node, join together children with ',' and add
            # outer parens
            child_strs = []
            cur = self.next
            while cur != self:
                child_strs.append(str(cur.branch))
                cur = cur.next
            tree_str = "(" + ",".join(child_strs) + ")" + branch_len_str

        if self.branch is None:
            # this is root, add semi-colon
            return tree_str + ";"

        return tree_str

        


def parse_branch_len(s):
    if not s.startswith(':'):
        raise ValueError("expected branch len string to start with ':'")

    # branch len can end with ',' or with ')'
    i1 = s.find(')')
    i2 = s.find(',')

    if (i1 == -1) or (i2 != -1 and i2 < i1):
        end = i2
    else:
        end = i1

    if end == -1:
        raise ValueError("expected branch length string '%s' to "
                         "contain a terminating ')' or ':'" % s)

    return float(s[1:end])
    


def parse_leaf_node(node_str, parent):
    node = TreeNode()
    node.branch = parent

    # leaf nodes are just connected to themselves
    node.next = node

    i = 0
    while (i < len(node_str)) and (node_str[i] not in (':', ',', ')')):
        i += 1

    if i == len(node_str):
        raise ValueError("expected name of leaf '%s' to have terminating ':'",
                         "',' or ')'" % node_str)

    node.name = node_str[0:i]

    if node_str[i] == ':':
        # there is a branch length, parse it too
        node.branch_len = parse_branch_len(node_str[i:])
    else:
        node.branch_len = None
                     
    return node
    

def parse_interior_node(node_str, parent):
    if len(node_str) < 2:
        raise ValueError("invalid node string")

    node = TreeNode()
    first_node = node
    node.branch = parent

    if not node_str.startswith('('):
        raise ValueError("node string '%s' should start with '('",
                         node_str)

    new_node_start = True

    end_idx = len(node_str)-1
    
    paren_count = 1
    for i in range(1, len(node_str)):
        if new_node_start:
            node.next = TreeNode()
            node = node.next

            if node_str[i] == '(':
                # this is the start of another interior node, call
                # this function recursively
                (node.branch, x) = parse_interior_node(node_str[i:], node)
            else:
                # this is the start of a leaf node
                node.branch = parse_leaf_node(node_str[i:], node)

            node.branch_len = node.branch.branch_len
            new_node_start = False
        
        if node_str[i] == ")":
            paren_count -= 1
            if paren_count == 0:
                # end of this internal node
                end_idx = i+1
                break
        elif node_str[i] == "(":
            paren_count += 1
        elif node_str[i] == ",":
            if paren_count == 1:
                # this is the start of another child node
                new_node_start = True

    # check paren count
    if paren_count != 0:
        raise ValueError("mismatched parens in tree: '%s'", node_str)


    # tie together nodes that represent interior node of this tree
    node.next = first_node

    if (end_idx < len(node_str)) and (node_str[end_idx] == ':'):
        # this internal node has a branch length
        first_node.branch_len = parse_branch_len(node_str[end_idx:])

    return (first_node, end_idx)
        


def parse_newick(newick_str):
    s = newick_str.strip()

    if s.endswith(";"):
        s = s[:-1]
    
    (node, end_idx) = parse_interior_node(s, None)

    if end_idx < len(s):
        raise ValueError("only parsed portion '%s' of complete tree "
                         "string '%s' missing paren or extra chars "
                         "(len:%d, end_idx:%d)?" %
                         (str(node), s, len(s), end_idx))
    
    return node





def main():
    sys.stderr.write("parsing tree:\n")
    tree = parse_newick("(((human:0.1,chimp:0.3):0.2,gorilla:0.1):0.1,(macaque:1.0,vervet:1.3):2.0)")
    sys.stderr.write("  %s\n" % str(tree))

    sys.stderr.write("copying tree:\n")
    tree = tree.copy()
    sys.stderr.write("  %s\n" % str(tree))    

    sys.stderr.write("finding node 'gorilla':\n")
    node = tree.find_node('gorilla')
    sys.stderr.write("  %s\n" % str(node))

    sys.stderr.write("removing gorilla:\n")
    tree = tree.remove_leaf('gorilla')
    sys.stderr.write("  %s\n" % str(tree))

    sys.stderr.write("removing chimp\n")
    tree = tree.remove_leaf('chimp')
    sys.stderr.write("  %s\n" % str(tree))

    sys.stderr.write("removing human\n")
    tree = tree.remove_leaf('human')
    sys.stderr.write("  %s\n" % str(tree))

    sys.stderr.write("\n\nparsing tree:\n")    
    tree = parse_newick("(((human:0.1,chimp:0.3):0.2,gorilla:0.1):0.1,(macaque:1.0,vervet:1.3):2.0)")
    sys.stderr.write("  %s\n" % str(tree))
    sys.stderr.write("removing human, vervet, chimp:\n")
    tree = tree.remove_leaves(['human', 'vervet', 'chimp'])
    sys.stderr.write("  %s\n" % str(tree))


if __name__ == "__main__":
    main()
    
