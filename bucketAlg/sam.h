#pragma once
#include <vector>
#include <unordered_map>
#include <cassert>
#include <string>
#include <stack>

class MySAM
{
public:
	struct State {
		int len, parent;
		int childPos;

		State() :len(0), parent(-1) {}
		State(int len) :len(len), parent(-1) {}

		std::unordered_map<int, int> next;
		
		std::vector<int> children;
	};
	std::vector<State> statePool;

	//at curState, given transition character, return the target state
	//if the transition does not exist, return -1
	inline int delta(int curState, int transition) {
		//assert( (curState, transition) in deltaMap )
		return statePool[curState].next[transition];
		//return deltaMap[std::make_pair(curState, transition)];
	}
	inline bool isDeltaExist(int curState, int trasition) {
		return statePool[curState].next.find(trasition) != statePool[curState].next.end();
	}
	inline void setDelta(int curState, int transition, int targetState) {
		statePool[curState].next[transition] = targetState;
		//deltaMap[std::make_pair(curState, transition)] = targetState;
	}
	inline int& len(int idx) {
		return statePool[idx].len;
	}
	inline int& parent(int idx) {
		return statePool[idx].parent;
	}
	inline std::unordered_map<int, int>& next(int idx) {
		return statePool[idx].next;
	}
	inline std::vector<int>& children(int idx) {
		return statePool[idx].children;
	}
	inline int& childPos(int idx) {
		return statePool[idx].childPos;
	}

	
	void setparent(int idx, int to) {
		//remove idx from the children of parent(idx)
		if (parent(idx) != -1) {
			childPos(children(parent(idx)).back()) = childPos(idx);
			std::swap(children(parent(idx))[childPos(idx)], children(parent(idx)).back());
			children(parent(idx)).pop_back();
		}
		//set parent(idx) 
		parent(idx) = to;
		children(to).push_back(idx);
		childPos(idx) = children(to).size() - 1;
	}

	int last;

	MySAM(int preAlloc=1) :last(0) {
		statePool.reserve(preAlloc);
		statePool.emplace_back(0);	//as root 
	}

	int extend(int c) {
		int cur = statePool.size();
		statePool.emplace_back(len(last) + 1);
	
		int p = last;
		//st[p].next.count(c) ->  delta(p, c) exists
		while (p != -1 && !isDeltaExist(p, c)) {
			setDelta(p, c, cur);	//let p--c-->cur
			p = parent(p);
		}
		if (p == -1) {
			//parent(cur) = 0;
			setparent(cur, 0);
		}
		else {
			int q = delta(p, c);
			if (len(p) + 1 == len(q)) {
				//parent(cur) = q;
				setparent(cur, q);
			}
			else {
				int clone = statePool.size();
				//statePool.emplace_back(len(p) + 1, parent(q));
				statePool.emplace_back(len(p) + 1);
				setparent(clone, parent(q));
				next(clone) = next(q);
				while (p != -1 && delta(p, c) == q) {
					setDelta(p, c, clone);
					p = parent(p);
				}
				//parent(cur) = clone;
				//parent(q) = clone;
				setparent(cur, clone);
				setparent(q, clone);
			}
		}
		last = cur;
		return last;
	}

	int extends(const std::initializer_list<int> &cs) {
		for (auto c : cs) {
			extend(c);
		}
		return last;
	}

	template<typename XS>
	int extends(const XS &cs) {
		for (auto c : cs) {
			extend(c);
		}
		return last;
	}

	//traverse antomaton, by default start from root
	//fmatch:: Idx->idx in xs->IO
	//fparent:: Idx->idx in xs->IO
	template<typename XS,typename FMATCH, typename Fparent>
	void traverse(const XS& xs, const FMATCH& fmatch, const Fparent& fparent, int curState=0) {
		int i = 0;
		for (auto& x : xs) {
			while (curState != 0 && !isDeltaExist(curState, x)) {
				curState = parent(curState);
				fparent(curState, i);
			}
			if (isDeltaExist(curState, x)) {
				//there exists such transition x
				curState = delta(curState, x);
				fmatch(curState, i);
			}
			++i;
		}
	}

	//f:: int -> bool (with some side effect probably)
	//apply f to all leaves of curState using dfs 
	//continue if f returns true
	//return if f returns false
	template<typename F>
	void mapSubtree(const F& f, int curState) {
		std::stack<int> stk;
		stk.push(curState);
		while (!stk.empty()) {
			int cur = stk.top();
			stk.pop();

			bool toContinue = f(cur);
			if (!toContinue) {
				return;
			}
			for (int c : children(cur)) {
				stk.push(c);
			}
		}
	}

	template<typename OS>
	void print(OS& os) {
		for (int i = 0; i < statePool.size();i++) {
			State& s = statePool[i];
			os << i << "(" << len(i) << ", [" << parent(i) << ", " << childPos(i) << "])" << "\n";
			os << "    children: ";
			for (int c : children(i)) {
				os << c << ", ";
			}
			os << "\n";
			for (auto&p : s.next) {
				os << i << "--" << p.first << "-->" << p.second << "\n";
			}
		}
	}
};