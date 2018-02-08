#ifndef EXPLORE_H_
#define EXPLORE_H_

template <
    typename State,
    typename Callback>
void explore(State & S, Callback callback = Callback()) {
  struct explorer {
    State & S;
    Callback callback;
    
    explorer(State & S, Callback const & callback)
        : S{S},
          callback{callback} {
    }
    
    bool explore() {
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
  
  explorer e{S, callback};
  e.explore();
}

#endif  // EXPLORE_H_
