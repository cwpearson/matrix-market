#pragma once

#include <algorithm>
#include <vector>

class KD {
public:
    struct Point {
        int i;
        int j;

        Point() = default;
        Point(int _i, int _j) : i(_i), j(_j) {}

        static bool by_ij(const Point &a, const Point &b) {
            if (a.i < b.i) {
                return true;
            } else if (a.i > b.i) {
                return false;
            } else {
                return a.j < b.j;
            }
        }

        static bool by_ji(const Point &a, const Point &b) {
            if (a.j < b.j) {
                return true;
            } else if (a.j > b.j) {
                return false;
            } else {
                return a.i < b.i;
            }
        }

    };
private:
    struct Node {
        Point location;
        Node *left;
        Node *right;
        Node() : left(nullptr), right(nullptr) {}
        ~Node() {
            if (left) delete left;
            if (right) delete right;
        }
    };

    Node *root_;

    Node *helper(Point *begin, Point *end, int depth) {
        // all work is done within the memory from begin to end

        if (begin >= end) {return nullptr;}

        const int iAxis = depth % 2;


        if (iAxis) {
            std::sort(begin, end, Point::by_ij);
        } else {
            std::sort(begin, end, Point::by_ji);
        }

        // split across median
        int mi = (end-begin) / 2;

        Node *n = new Node;

        // split points across median
        n->location = begin[mi];
        n->left = helper(begin, begin+mi, depth+1);
        n->right = helper(begin+mi+1, end, depth+1);
        
        return n;
    }



    // count all points in [ilb...iub) and [jlb...jub)
    int range_count_helper(Node *n, int ilb, int iub, int jlb ,int jub, int depth) const {

        int i = n->location.i;
        int j = n->location.j;

        int count = 0;
        // node location is in the range
        if (i >= ilb && i < iub && j >= jlb && j < jub) {
            count += 1;
        }


        int iAxis = depth % 2;

        if (iAxis) {
            // TODO: optimization if subtree is totally contained
            if (i >= ilb && n->left) {
                count += range_count_helper(n->left, ilb, iub, jlb, jub, depth+1);
            }
            if (i < iub && n->right) {
                count += range_count_helper(n->right, ilb, iub, jlb, jub, depth+1);
            }
        } else {
            // TODO: optimization if subtree is totally contained
            if (j >= jlb && n->left) {
                count += range_count_helper(n->left, ilb, iub, jlb, jub, depth+1);
            }
            if (j < jub && n->right) {
                count += range_count_helper(n->right, ilb, iub, jlb, jub, depth+1);
            }
        }
        return count;
    }

public:
    KD(std::vector<Point> ps /* by value so we can use it as scratch*/) {
        root_ = helper(&ps[0], &ps[ps.size()], 0);
    }
    ~KD() {
        delete root_;
    }

    // find all points in [ilb...iub) and [jlb...jub)
    int range_count(int ilb, int iub, int jlb ,int jub) const {
        return range_count_helper(root_, ilb, iub, jlb , jub, 0);
    }

};