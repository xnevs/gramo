#ifndef EXPLORE_H_
#define EXPLORE_H_

template <
    typename State,
    typename Callback>
void explore(State & S, Callback callback = Callback()) {
  int count = 0;
  struct explorer {
    State & S;
    Callback callback;
    
    int & count;
    
    explorer(State & S, Callback const & callback, int & count)
        : S{S},
          callback{callback},
          count{count} {
    }
    
    bool explore() {
      ++count;
      if (S.full()) {
        return callback(S);
      } else {
        bool proceed = true;
        for (auto y : S.candidates()) {
          bool proceed = true;
          S.advance();
          bool success = S.assign(y);
          if (success) {
            S.push(y);
            proceed = explore();
            S.pop();
          }
          S.revert();
          if (!proceed) {
            break;
          }
        }
        return proceed;
      }
    }
  };
  
  explorer e{S, callback, count};
  e.explore();
  std::cout << "count: " << count << std::endl;
}

#endif  // EXPLORE_H_
