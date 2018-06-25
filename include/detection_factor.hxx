#ifndef LP_MP_DETECTION_FACTOR_HXX
#define LP_MP_DETECTION_FACTOR_HXX

#include "config.hxx"
#include <bitset>
#include "vector.hxx"
#include "sat_solver.hxx"

namespace LP_MP {

enum class exit_constraint_position {lower,upper}; // detection is non-overlapping below/above some other

// factor containing whether cell hypothesis is true, all its possible predecessors and successors.

// primal has the following meaning:
// incoming_edge_ = std::numeric_limits<INDEX>::max() means no decision is taken yet.
//                = std::numeric_limits<REAL>::max()-1 means that no primal is to be taken.
//                < incoming.size() means incoming edge is assigned
// same for outgoing_edge_
class detection_factor {

  friend class transition_message;
  friend class at_most_one_cell_message;
  template<exit_constraint_position POSITION> friend class exit_constraint_message;

public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_edge_taken = std::numeric_limits<INDEX>::max()-1;

  detection_factor(
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost)
    : 
      detection(detection_cost),
      incoming(no_incoming_transition_edges + no_incoming_division_edges + 1, 0.0),
      outgoing(no_outgoing_transition_edges + no_outgoing_division_edges + 1, 0.0),
      min_incoming_dirty_(true),
      min_outgoing_dirty_(true)
  {
      std::fill(incoming.begin(), incoming.end(), std::numeric_limits<REAL>::infinity());
      std::fill(outgoing.begin(), outgoing.end(), std::numeric_limits<REAL>::infinity());
    incoming[ no_incoming_transition_edges + no_incoming_division_edges ] = appearance_cost;
    outgoing[ no_outgoing_transition_edges + no_outgoing_division_edges ] = disappearance_cost;
    //std::cout << "no incoming = " << incoming.size() << ", no outgoing = " << outgoing.size() << "\n";
  }

  REAL min_detection_cost() const 
  {
    //std::cout << pot_[0] << " : ";
    //for(auto it=incoming.begin(); it!= incoming.end(); ++it) std::cout << *it << ",";
    //std:cout << " : ";
    //for(auto it=outgoing.begin(); it!= outgoing.end(); ++it) std::cout << *it << ",";
    //std::cout << " = ";

    assert(incoming.size() > 0);
    assert(outgoing.size() > 0);
    return detection + min_incoming() + min_outgoing();
  }
  
  void MaximizePotentialAndComputePrimal() 
  {
    if(incoming_edge_ < incoming.size() && outgoing_edge_ < outgoing.size()) {
      return;
    }
    if(incoming_edge_ == no_edge_taken && incoming.size() > 0) {
      outgoing_edge_ = no_edge_taken;
      return;
    }
    if(outgoing_edge_ == no_edge_taken && outgoing.size() > 0) {
      incoming_edge_ = no_edge_taken;
      return;
    }

    INDEX incoming_cand = no_edge_taken;
    REAL lb = detection;
    if(incoming.size() > 0) {
      incoming_cand = std::min_element(incoming.begin(), incoming.end()) - incoming.begin();
      lb += incoming[incoming_cand]; 
    }
    INDEX outgoing_cand = no_edge_taken;
    if(outgoing.size() > 0) {
      outgoing_cand = std::min_element(outgoing.begin(), outgoing.end()) - outgoing.begin();
      lb += outgoing[outgoing_cand]; 
    }
    // one edge already labelled
    if(incoming_edge_ < incoming.size() && outgoing.size() > 0 && outgoing_edge_ == no_primal_decision) {
      assert(outgoing_edge_ == no_primal_decision);
      outgoing_edge_ = outgoing_cand;
      return;
    } else if(outgoing_edge_ < outgoing.size() && incoming.size() > 0 && incoming_edge_ == no_primal_decision) {
      assert(incoming_edge_ == no_primal_decision);
      incoming_edge_ = incoming_cand;
      return;
    }
    // no edge labelled yet
    if(lb < 0.0) {
      incoming_edge_ = incoming_cand;
      outgoing_edge_ = outgoing_cand;
      return;
    } else {
      incoming_edge_ = no_edge_taken;
      outgoing_edge_ = no_edge_taken;
    }
  }
  REAL LowerBound() const {
      // check whether every entry has been set
      //std::cout << detection << "; ";
      //for(auto x : incoming) { std::cout << x << ","; }
      //std::cout << ";";
      //for(auto x : outgoing) { std::cout << x << ","; }
      //std::cout << "\n";
      assert(*std::max_element(incoming.begin(), incoming.end()) < std::numeric_limits<REAL>::infinity());
      assert(*std::max_element(outgoing.begin(), outgoing.end()) < std::numeric_limits<REAL>::infinity());
      return std::min(min_detection_cost(), REAL(0.0));
  }

  //REAL& detection_cost() { return *pot_; }

  /*
  REAL* incoming.begin() const { return incoming_.begin(); };
  REAL* incoming.end() const { return incoming_.end(); };
  INDEX incoming.size() const { return incoming_.size(); }

  REAL* outgoing.begin() const { return outgoing_.begin(); };
  REAL* outgoing.end() const { return outgoing_.end(); };
  INDEX outgoing.size() const { return outgoing_.size(); }
  */

  REAL EvaluatePrimal() const {
    assert(incoming.size() > 0);
    assert(outgoing.size() > 0);
    if(incoming_edge_ < incoming.size() && outgoing_edge_ < outgoing.size()) {
      return detection + incoming[incoming_edge_] + outgoing[outgoing_edge_];
    } else {
      return 0.0;
    }
  }

  void set_incoming_transition_cost(const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(incoming[edge_index] == std::numeric_limits<REAL>::infinity()); 
    incoming[edge_index] = cost; 
  } 
  void set_outgoing_transition_cost(const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(outgoing[edge_index] == std::numeric_limits<REAL>::infinity());
    outgoing[edge_index] = cost;
  }
  void set_incoming_division_cost(const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(incoming[no_incoming_transition_edges + edge_index] == std::numeric_limits<REAL>::infinity());
    incoming[no_incoming_transition_edges + edge_index] = cost; 
  } 
  void set_outgoing_division_cost(const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(outgoing[no_outgoing_transition_edges + edge_index] == std::numeric_limits<REAL>::infinity());
    outgoing[no_outgoing_transition_edges + edge_index] = cost; 
  }

  REAL detection_cost() const { return detection; }
  REAL appearance_cost() const { return incoming[incoming.size()-1]; }
  REAL disappearance_cost() const { return outgoing[outgoing.size()-1]; }

  INDEX size() const { return 1 + incoming.size() + outgoing.size(); }

  void init_primal() 
  { 
    incoming_edge_ = no_primal_decision;
    outgoing_edge_ = no_primal_decision; 
  }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( incoming_edge_, outgoing_edge_ ); }
  //template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cereal::binary_data( pot_, sizeof(REAL)*size()) ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( detection, incoming, outgoing ); }

  auto export_variables() { return std::tie( detection, incoming, outgoing ); }

  INDEX outgoing_edge() const { return outgoing_edge_; }
  INDEX& outgoing_edge() { return outgoing_edge_; }
  INDEX incoming_edge() const { return incoming_edge_; }
  INDEX& incoming_edge() { return incoming_edge_; }

  bool detection_active() const {
    return (incoming_edge_ < incoming.size() || outgoing_edge_ < outgoing.size());
  }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::variable detection_var, typename SOLVER::vector incoming_var, typename SOLVER::vector outgoing_var) const
  {
    // create variables
    auto outgoing_sum = s.add_at_most_one_constraint(outgoing_var.begin(), outgoing_var.end());
    auto incoming_sum = s.add_at_most_one_constraint(incoming_var.begin(), incoming_var.end());

    // detection var must be equal to incoming and outgoing var
    s.make_equal(detection_var,incoming_sum);
    s.make_equal(detection_var,outgoing_sum);
  }

  /*
  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    sat_literal detection_literal;
    sat_literal_vector incoming_literals(incoming.size());
    sat_literal_vector outgoing_literals(outgoing.size());
    load_sat_literals(begin, detection_literal, incoming_literals, outgoing_literals);

    const REAL cost = min_detection_cost();
    if(cost <= -th) {
      //std::cout << "in reduction: detection factor must be on\n";
      assumptions.push_back(detection_literal);
    }
    if(cost <= th) {
      const REAL incoming_min = min_incoming(); //*std::min_element(incoming.begin(), incoming.end());
      for(INDEX i=0; i<incoming.size(); ++i) {
        if(incoming[i] > incoming_min + th) { 
           assumptions.push_back(-incoming_literals[i]);
         }
      }

      const REAL outgoing_min = min_outgoing(); //*std::min_element(outgoing.begin(), outgoing.end());
      for(INDEX i=0; i<outgoing.size(); ++i) {
        if(outgoing[i] > outgoing_min + th) { 
           assumptions.push_back(-outgoing_literals[i]);
         }
      } 
    } else {
      assumptions.push_back(-detection_literal);
    }
  }
  */

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::variable detection_var, typename SOLVER::vector incoming_var, typename SOLVER::vector outgoing_var)
  {
    if(s.solution(detection_var)) {
      incoming_edge_ = s.first_active(incoming_var);
      outgoing_edge_ = s.first_active(outgoing_var);
    } else {
      incoming_edge_ = no_edge_taken;
      outgoing_edge_ = no_edge_taken;
    }
  }

  // invoke in general case
  void update_incoming(const INDEX edge_index, const REAL delta)
  {
    const auto prev_val = incoming[edge_index];
    incoming[edge_index] += delta;
    if(min_incoming_dirty_) {
      return; 
    } else {
      if(prev_val == min_incoming_ && delta > 0) {
        min_incoming_dirty_ = true;
        return;
      } else {
        min_incoming_ = std::min(min_incoming_, incoming[edge_index]);
        return;
      } 
    }
  }

  // invoke when we have have computed delta in terms of this factor. We do not need to invalidate min_incoming in any case
  void update_incoming_valid(const INDEX edge_index, const REAL delta)
  {
    return update_incoming(edge_index, delta);
    //if(incoming[edge_index] == min_incoming_) { min_incoming_ += delta; }
    incoming[edge_index] += delta; 
    //if(!min_incoming_dirty_) { assert(min_incoming_ == incoming.min()); }
  }

  // invoke in general case
  void update_outgoing(const INDEX edge_index, const REAL delta)
  {
    const auto prev_val = outgoing[edge_index];
    outgoing[edge_index] += delta;
    if(min_outgoing_dirty_) {
      return; 
    } else {
      if(prev_val == min_outgoing_ && delta > 0) {
        min_outgoing_dirty_ = true;
        return;
      } else {
        min_outgoing_ = std::min(min_outgoing_, outgoing[edge_index]);
        return;
      } 
    }
  }

  // invoke when we have have computed delta in terms of this factor. We do not need to invalidate min_outgoing in any case
  void update_outgoing_valid(const INDEX edge_index, const REAL delta)
  {
    return update_outgoing(edge_index, delta);
    //if(outgoing[edge_index] == min_outgoing_) { min_outgoing_ += delta; }
    outgoing[edge_index] += delta; 
    //if(!min_outgoing_dirty_) { assert(min_outgoing_ == outgoing.min()); }
  }


  REAL min_incoming() const
  {
    //return incoming.min();
    if(min_incoming_dirty_) {
      min_incoming_ = incoming.min();
      min_incoming_dirty_ = false;
    }
    assert(min_incoming_ == incoming.min());
    return min_incoming_;
  }

  REAL min_outgoing() const
  {
    //return outgoing.min();
    if(min_outgoing_dirty_) {
      min_outgoing_ = outgoing.min();
      min_outgoing_dirty_ = false;
    }
    assert(min_outgoing_ == outgoing.min());
    return min_outgoing_;
  } 

  REAL min_incoming_marginal_diff(const INDEX edge_index) const
  {
    const REAL detection_outgoing_cost = detection + min_outgoing();

    const REAL incoming_min = min_incoming();
    const REAL incoming_val = incoming[edge_index];
    assert(incoming_val >= incoming_min);

    if(incoming_val != incoming_min) {
      return detection_outgoing_cost + incoming_val - std::min(detection_outgoing_cost + incoming_min, REAL(0.0)); 
    } else {
      incoming[edge_index] = std::numeric_limits<REAL>::infinity();
      const REAL second_incoming_min = incoming.min();
      incoming[edge_index] = incoming_val;

      return detection_outgoing_cost + incoming_val - std::min(detection_outgoing_cost + second_incoming_min, REAL(0.0)); 
    }
  } 

  REAL min_outgoing_marginal_diff(const INDEX edge_index) const
  {
    const REAL detection_incoming_cost = detection + min_incoming();

    const REAL outgoing_min = min_outgoing();
    const REAL outgoing_val = outgoing[edge_index];
    assert(outgoing_val >= outgoing_min);

    if(outgoing_val != outgoing_min) {
      return detection_incoming_cost + outgoing_val - std::min(detection_incoming_cost + outgoing_min, REAL(0.0));
    } else {
      outgoing[edge_index] = std::numeric_limits<REAL>::infinity();
      const REAL second_outgoing_min = outgoing.min();
      outgoing[edge_index] = outgoing_val;

      return detection_incoming_cost + outgoing_val - std::min(detection_incoming_cost + second_outgoing_min, REAL(0.0));
    }
  }

  REAL detection;
  //mutable small_vector<REAL,16> incoming;
  //mutable small_vector<REAL,16> outgoing;
  mutable vector<REAL> incoming;
  mutable vector<REAL> outgoing;

private:
  // for efficient computation of messages, keep track of smallest and second smallest element in incoming and outgoing vector
  mutable REAL min_incoming_, min_outgoing_;

  // primal  
  INDEX incoming_edge_, outgoing_edge_;

  mutable bool min_incoming_dirty_ = true, min_outgoing_dirty_ = true; 
};

class mapping_edge_factor {
public:
  mapping_edge_factor(const REAL c) : cost(c) {}

  REAL LowerBound() const
  { 
    return std::min(REAL(0.0), cost);
  }

  REAL EvaluatePrimal() const
  {
    if(primal_) {
      return cost;
    } else {
      return 0.0;
    }
  }
       
  INDEX size() const { return 1; }

  void init_primal() {}
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cost ); }
  auto export_variables() { return std::tie(cost); }

  REAL cost;

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::variable) const
  {}

  /*
  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    if(cost <= -th) {
      assumptions.push_back(begin);
    } else if(cost >= th) {
      assumptions.push_back(-begin);
    }
  }
  */

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::variable var)
  {
    //assert(s.get_model()[first] != CMSat::l_Undef);
    //std::cout << lglmaxvar(s) << "\n";
    primal_ = s.solution(var);
  }

private: 
  bool primal_;
};

// left is mapping_edge_factor, right is detection_factor
class cell_incoming_edge_detection_factor {
public:
  cell_incoming_edge_detection_factor(const INDEX incoming_edge_index) : incoming_edge_index_(incoming_edge_index) {}

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    l.cost += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 1 || msg_dim == 2);
    if(msg_dim == 1) {
      r.update_incoming(incoming_edge_index_, msg);
    } else {
      r.update_incoming_valid(incoming_edge_index_, msg);
    }
    
    //assert(msg_dim == 0);
    //r.incoming[incoming_edge_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg[1] -= omega*l.cost;
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    msg[2] -= omega*r.min_incoming_marginal_diff(incoming_edge_index_);
  } 

  // send messages from detection factor along outgoing edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    //std::cout << "send incoming with omega = " << omega << "\n";

    const REAL detection_outgoing_cost = r.detection + r.min_outgoing();

    std::vector<bool> edge_taken(r.incoming.size(), false);
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
      edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = true;
    }

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    const auto incoming_size = r.incoming.size();
    for(INDEX i=0; i<incoming_size; ++i) {
      const REAL val = r.incoming[i];
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    const REAL set_to_cost = std::min(detection_outgoing_cost + std::min(second_smallest_taken, smallest_not_taken), REAL(0.0));

    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const REAL msg = detection_outgoing_cost + r.incoming[ (*msg_it).GetMessageOp().incoming_edge_index_ ] - set_to_cost;
      (*msg_it)[0] -= omega*msg;
    } 
    return;

    /*
    const auto smallest_incoming = two_smallest_elements<REAL>(r.incoming.begin(), r.incoming.end()); // do not take into account disappearance cost 

    //assert(false); // wrong message computation?
    const REAL set_to_cost = std::min(detection_outgoing_cost + smallest_incoming[1], REAL(0.0));
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const REAL msg = detection_outgoing_cost + r.incoming[ (*msg_it).GetMessageOp().incoming_edge_index_ ] - set_to_cost;
      (*msg_it)[1] -= omega*msg;
    } 
    */
  }


  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_var,
      RIGHT_FACTOR& r, typename SOLVER::variable detection_var, typename SOLVER::vector incoming_vars, typename SOLVER::vector outgoing_vars) const
  {
    s.make_equal(left_var, incoming_vars[incoming_edge_index_]);
  }

private:
  const INDEX incoming_edge_index_; 
};

// left is mapping_edge_factor, right is detection_factor
class cell_outgoing_edge_detection_factor {
public:
  cell_outgoing_edge_detection_factor(const INDEX outgoing_edge_index, const bool split)
    : outgoing_edge_index_(outgoing_edge_index),
    split_(split)
  {}

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    //assert(msg_dim == 0);
    l.cost += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 1 || msg_dim == 2);
    if(msg_dim == 1) {
      r.update_outgoing_valid(outgoing_edge_index_, msg);
    } else {
      r.update_outgoing(outgoing_edge_index_, msg);
    }

    //assert(msg_dim == 0);
    //r.outgoing[outgoing_edge_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg[2] -= omega*l.cost;
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    msg[1] -= omega*r.min_outgoing_marginal_diff(outgoing_edge_index_);
  } 

  // send messages from detection factor along incoming edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    //std::cout << "send outgoing with omega = " << omega << "\n";

    const REAL detection_incoming_cost = r.detection + r.min_incoming();

    std::vector<bool> edge_taken(r.outgoing.size(), false);
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
      edge_taken[(*msg_it).GetMessageOp().outgoing_edge_index_] = true;
    }

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    const auto outgoing_size = r.outgoing.size();
    for(INDEX i=0; i<outgoing_size; ++i) {
      const REAL val = r.outgoing[i];
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    const REAL set_to_cost = std::min(detection_incoming_cost + std::min(second_smallest_taken, smallest_not_taken), REAL(0.0));
    for(; msg_begin!=msg_end; ++msg_begin) {
        const REAL w = (*msg_begin).GetMessageOp().split_ ? 0.5 : 1.0;
        const REAL msg = w*(detection_incoming_cost + r.outgoing[ (*msg_begin).GetMessageOp().outgoing_edge_index_ ] - set_to_cost);
        (*msg_begin)[0] -= omega*msg;
    } 

    return;

    /*
    // compute smallest and second smallest value over all outgoing edges
    const auto smallest_outgoing = two_smallest_elements<REAL>(r.outgoing.begin(), r.outgoing.end());

    //assert(false); // wrong message computation?
    const REAL set_to_cost = std::min(detection_incoming_cost + smallest_outgoing[1], REAL(0.0));
    for(; msg_begin!=msg_end; ++msg_begin) {
      const REAL w = (*msg_begin).GetMessageOp().split_ ? 0.5 : 1.0;
      const REAL msg = w*(detection_incoming_cost + r.outgoing[ (*msg_begin).GetMessageOp().outgoing_edge_index_ ] - set_to_cost);
      (*msg_begin)[1] -= omega*msg;
    } 
    return; 
    */
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_var,
      RIGHT_FACTOR& r, typename SOLVER::variable detection_var, typename SOLVER::vector incoming_vars, typename SOLVER::vector outgoing_vars) const
  {
    s.make_equal(left_var, outgoing_vars[outgoing_edge_index_]);
  }

private:
  // to do: use bitfields
  const INDEX outgoing_edge_index_; 
  const bool split_;
}; 

// message connecting outgoing edge to incoming edge of detection factors between consecutive timeframes
class transition_message {
public:

  transition_message(const bool split, const INDEX outgoing_edge_index, const INDEX incoming_edge_index) 
    : 
    outgoing_edge_index_(outgoing_edge_index),
    incoming_edge_index_(incoming_edge_index),
    split_(split)
  {
    assert(outgoing_edge_index <= std::numeric_limits<SHORT_INDEX>::max());
    assert(incoming_edge_index <= std::numeric_limits<SHORT_INDEX>::max());
  }

  // send messages from detection factor along outgoing edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    //std::cout << "send transition messages to left with omega = " << omega << "\n";

    // check #messages+1 <= no incoming edges
#ifndef NDEBUG
    {
      INDEX c=0;
      for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
      assert(c+1 <= rightFactor.incoming.size());
    }
#endif 

    const REAL detection_outgoing_cost = rightFactor.detection + rightFactor.min_outgoing();

    //std::vector<bool> edge_taken(rightFactor.incoming.size(), false);
    //std::bitset<128> edge_taken(false);
    std::array<unsigned char,128> edge_taken;
    std::fill(edge_taken.begin(), edge_taken.end(), 0);
    assert(rightFactor.incoming.size() <= edge_taken.size());
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
        const auto idx = (*msg_it).GetMessageOp().incoming_edge_index_;
        assert(edge_taken[idx] == 0);
        edge_taken[idx] = 1;
    }

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    const auto incoming_size = rightFactor.incoming.size();
    for(INDEX i=0; i<incoming_size; ++i) {
      const auto val = rightFactor.incoming[i];
      if(edge_taken[i]) {
        const auto min = std::min(smallest_taken, val);
        const auto max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }
    assert(std::min(smallest_taken, smallest_not_taken) == rightFactor.min_incoming());

    const auto set_to_cost = std::min(detection_outgoing_cost + std::min(smallest_not_taken, second_smallest_taken), REAL(0.0));

    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
        const auto idx = (*msg_it).GetMessageOp().incoming_edge_index_;
        const auto msg = detection_outgoing_cost + rightFactor.incoming[ idx ] - set_to_cost;
        (*msg_it)[0] -= omega*msg;
    } 
  }

  template<typename RIGHT_FACTOR, typename MSG_ITERATOR>
  REAL send_messages_to_left_improvement(RIGHT_FACTOR& r, MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
  {
    const REAL detection_outgoing_cost = r.detection + r.min_outgoing();

    //std::vector<bool> edge_taken(rightFactor.incoming.size(), false);
    //std::bitset<128> edge_taken(false);
    std::array<unsigned char,128> edge_taken;
    std::fill(edge_taken.begin(), edge_taken.end(), 0);
    assert(r.incoming.size() <= edge_taken.size());
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
      edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = 1;
    }

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    const auto incoming_size = r.incoming.size();
    for(INDEX i=0; i<incoming_size; ++i) {
      const auto val = r.incoming[i];
      if(edge_taken[i]) {
        const auto min = std::min(smallest_taken, val);
        const auto max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    const auto right_lb = std::min(detection_outgoing_cost + r.min_incoming(), 0.0);
    const auto new_right_lb = std::min(detection_outgoing_cost + std::min(smallest_not_taken, second_smallest_taken), 0.0);

    REAL left_lb = 0.0;
    REAL new_left_lb = 0.0;
    const auto set_to_cost = std::min(detection_outgoing_cost + std::min(smallest_not_taken, second_smallest_taken), REAL(0.0));

    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const auto diff = detection_outgoing_cost + r.incoming[ (*msg_it).GetMessageOp().incoming_edge_index_ ] - set_to_cost;
      auto& l = *(*msg_it).GetLeftFactor()->GetFactor();
      left_lb += l.LowerBound();
      const auto outgoing_edge_index = (*msg_it).GetMessageOp().outgoing_edge_index_;
      const auto prev_val = l.outgoing[outgoing_edge_index];
      l.outgoing[outgoing_edge_index] += diff;
      const REAL new_outgoing_cost = l.outgoing.min();
      l.outgoing[outgoing_edge_index] = prev_val;
      new_left_lb += std::min(0.0, l.min_incoming() + l.detection + new_outgoing_cost); 
    } 

    assert(new_left_lb + new_right_lb - (left_lb + right_lb) >= -eps); 
    return new_left_lb + new_right_lb - (left_lb + right_lb); 
  }

  // send messages from detection factor along incoming edges
  template<typename LEFT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToRight(const LEFT_FACTOR& leftFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    //std::cout << "send transition messages to right with omega = " << omega << "\n";

    const REAL detection_incoming_cost = leftFactor.detection + leftFactor.min_incoming();

    //std::vector<bool> edge_taken(leftFactor.outgoing.size(), false);
    //std::bitset<128> edge_taken(false);
    std::array<unsigned char,128> edge_taken;
    assert(leftFactor.outgoing.size() <= edge_taken.size());
    std::fill(edge_taken.begin(), edge_taken.end(), 0);
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
        const auto idx = (*msg_it).GetMessageOp().outgoing_edge_index_;
        edge_taken[idx]++;
    }

    vector<REAL> outgoing_copy(leftFactor.outgoing); // we need an outgoing copy, since for shared outgoing edges (i.e. branching) we need the value of outgoing twice, yet it may be changed by reparametrization

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    const auto outgoing_size = outgoing_copy.size();
    auto val_it = outgoing_copy.begin();
    for(std::size_t i=0; i<outgoing_size; ++i, ++val_it) {
      const REAL val = *val_it;
      if(edge_taken[i] > 0) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }
    assert(std::min(smallest_taken, smallest_not_taken) == leftFactor.min_outgoing());

    const REAL set_to_cost = std::min(detection_incoming_cost + std::min(smallest_not_taken, second_smallest_taken), REAL(0.0));
     
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const auto idx = (*msg_it).GetMessageOp().outgoing_edge_index_ ;
      assert(edge_taken[idx] > 0);
      const REAL w = 1.0/REAL(edge_taken[idx]);
      //const REAL w = (*msg_it).GetMessageOp().split_ ? 0.5 : 1.0;
      const REAL msg = w*(detection_incoming_cost + outgoing_copy[ idx ] - set_to_cost);
      (*msg_it)[0] -= omega*msg;
    } 
  }

  template<typename LEFT_FACTOR, typename MSG_ITERATOR>
  REAL send_messages_to_right_improvement(LEFT_FACTOR& l, MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
  {
    const REAL detection_incoming_cost = l.detection + l.min_incoming();

    std::array<unsigned char,128> edge_taken;
    std::fill(edge_taken.begin(), edge_taken.end(), 0);
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
      edge_taken[(*msg_it).GetMessageOp().outgoing_edge_index_] = 1;
    }
    assert(l.outgoing.size() <= edge_taken.size());

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    const auto outgoing_size = l.outgoing.size();
    for(INDEX i=0; i<outgoing_size; ++i) {
      const REAL val = l.outgoing[i];
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }


    const auto left_lb = std::min(detection_incoming_cost + l.min_outgoing(), 0.0);
    const auto new_left_lb = std::min(detection_incoming_cost + std::min(smallest_not_taken, second_smallest_taken), REAL(0.0));

    REAL right_lb = 0.0;
    REAL new_right_lb = 0.0;
    const REAL set_to_cost = std::min(detection_incoming_cost + std::min(smallest_not_taken, second_smallest_taken), REAL(0.0));
     
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const REAL w = (*msg_it).GetMessageOp().split_ ? 0.5 : 1.0;
      const auto diff = w*(detection_incoming_cost + l.outgoing[ (*msg_it).GetMessageOp().outgoing_edge_index_ ] - set_to_cost);
      auto& r = *(*msg_it).GetRightFactor()->GetFactor();
      right_lb += r.LowerBound();

      const auto incoming_edge_index = (*msg_it).GetMessageOp().incoming_edge_index_;
      // only do this for positive diff
      const auto prev_val = r.incoming[incoming_edge_index];
      r.incoming[incoming_edge_index] += diff;
      const REAL new_incoming_cost = r.incoming.min();
      r.incoming[incoming_edge_index] = prev_val; 
      new_right_lb += std::min(0.0, r.min_outgoing() + r.detection + new_incoming_cost);
    } 
    assert(new_left_lb + new_right_lb - (left_lb + right_lb) >= -eps); 
    return new_left_lb + new_right_lb - (left_lb + right_lb); 
  }

/*
  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    if(l.outgoing_edge_ == outgoing_edge_index_ && r.incoming_edge_ != incoming_edge_index_) {
      r.incoming_edge_ = incoming_edge_index_;
      return true;
    } else {
      return false;
    }
  }
  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) 
  {
    if(r.incoming_edge_ == incoming_edge_index_ && r.outgoing_edge_ != outgoing_edge_index_) {
      l.outgoing_edge_ = outgoing_edge_index_;
      return true;
    } else {
      return false;
    }
  }
  */

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    msg[0] -= omega*r.min_incoming_marginal_diff(incoming_edge_index_);
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  REAL send_message_to_left_improvement(LEFT_FACTOR& l, RIGHT_FACTOR& r) const
  {
      const auto right_lb = r.LowerBound();
      const auto diff = r.min_incoming_marginal_diff(incoming_edge_index_);
      const auto new_right_lb = std::min(0.0, r.min_incoming() + r.detection_cost() + std::min(r.min_outgoing(), r.outgoing[outgoing_edge_index_] + diff));
      if(diff >= 0) {
          return new_right_lb - right_lb;
      } else {
          const auto left_lb = l.LowerBound();
          const auto new_left_lb = std::min(0.0, r.min_outgoing() + r.detection_cost() + std::min(r.min_incoming(), r.incoming[incoming_edge_index_] - diff));
          return new_left_lb + new_right_lb - (left_lb + right_lb);
      }
  }

  template<typename RIGHT_FACTOR, typename G2>
  void ReceiveRestrictedMessageFromRight(RIGHT_FACTOR& r, G2& msg)
  { 
    if(r.incoming_edge_ < r.incoming.size()) { // incoming edge has already been picked
      if(r.incoming_edge_ == incoming_edge_index_) {
        const REAL val_prev = msg.GetLeftFactor()->GetFactor()->outgoing[outgoing_edge_index_];
        msg[0] -= -std::numeric_limits<REAL>::infinity(); 
        const REAL val = msg.GetLeftFactor()->GetFactor()->outgoing[outgoing_edge_index_];
        //assert(msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_) < -10000);
      } else {
        const REAL val_prev = msg.GetLeftFactor()->GetFactor()->outgoing[outgoing_edge_index_];
        msg[0] -= std::numeric_limits<REAL>::infinity(); 
        const REAL val = msg.GetLeftFactor()->GetFactor()->outgoing[outgoing_edge_index_];
        //assert(msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_) > 10000);
      }
    } else { // no incoming edge has been picked yet.
      return;
      //send_message_to_left(r,msg); 
    }
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg[0] -= omega* l.min_outgoing_marginal_diff(outgoing_edge_index_);
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  REAL send_message_to_right_improvement(LEFT_FACTOR& l, RIGHT_FACTOR& r) const
  {
      const auto left_lb = l.LowerBound();
      const auto diff = l.min_outgoing_marginal_diff(outgoing_edge_index_);
      const auto new_left_lb = std::min(0.0, r.min_outgoing() + r.detection_cost() + std::min(r.min_incoming(), r.incoming[incoming_edge_index_] - diff));
      if(diff >= 0) {
          return new_left_lb - left_lb;
      } else {
          const auto right_lb = l.LowerBound();
          const auto new_right_lb = std::min(0.0, r.min_incoming() + r.detection_cost() + std::min(r.min_outgoing(), r.outgoing[outgoing_edge_index_] + diff));
          return new_left_lb + new_right_lb - (left_lb + right_lb);
      }
  }

  template<typename LEFT_FACTOR, typename G2>
  void ReceiveRestrictedMessageFromLeft(LEFT_FACTOR& l, G2& msg)
  { 
    if(l.outgoing_edge_ < l.outgoing.size()) { // outgoing edge has already been picked
      if(l.outgoing_edge_ == outgoing_edge_index_) {
        msg[0] -= -std::numeric_limits<REAL>::infinity(); 
        assert(msg.GetRightFactor()->GetFactor()->incoming[incoming_edge_index_] < -10000);
      } else {
        msg[0] -= std::numeric_limits<REAL>::infinity(); 
        assert(msg.GetRightFactor()->GetFactor()->incoming[incoming_edge_index_] > 10000);
      }
    } else { // no incoming edge has been picked yet.
      return;
      //send_message_to_right(l,msg); 
    }
  }

  // msg dim = 1: update left valid, update right general
  // msg dim = 2: update right valid, update left general
  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
      assert(msg_dim == 0);
      l.update_outgoing(outgoing_edge_index_, msg);
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
      assert(msg_dim == 0);
      r.update_incoming(incoming_edge_index_, msg);
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    if(l.outgoing_edge_ == outgoing_edge_index_) {
      return r.incoming_edge_ == incoming_edge_index_;
    }
    if(r.incoming_edge_ == incoming_edge_index_) {
      return r.outgoing_edge_ == outgoing_edge_index_;
    } 
    return true;
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_detection_var, typename SOLVER::vector left_incoming_vars, typename SOLVER::vector left_outgoing_vars,
      RIGHT_FACTOR& r, typename SOLVER::variable right_detection_var, typename SOLVER::vector right_incoming_vars, typename SOLVER::vector right_outgoing_vars) const
  {
    assert(outgoing_edge_index_ < l.outgoing.size());
    assert(incoming_edge_index_ < r.incoming.size());
    s.make_equal(left_outgoing_vars[outgoing_edge_index_], right_incoming_vars[incoming_edge_index_]);
  }
private:
  const SHORT_INDEX outgoing_edge_index_;
  const SHORT_INDEX incoming_edge_index_;
  const bool split_; // is split really needed? Whether it is a split transition can be found out in SendMessagesTo... by checking whether there is two outgoing_edge_indices
};

// multiple cell detection hypotheses can be mutually exclusive. 
// simplex x1 + ... + xn <= 1
// to account for overlapping detections: only one can be active
class at_most_one_cell_factor : public vector<REAL> {

  friend class at_most_one_cell_message;

public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_primal_active = std::numeric_limits<INDEX>::max()-1;
  constexpr static INDEX primal_infeasible = std::numeric_limits<INDEX>::max()-2;

   at_most_one_cell_factor(const INDEX size) : vector<REAL>(size, 0.0)
   //at_most_one_cell_factor(const INDEX size) : small_vector<REAL,4>(size, 0.0)
   {}

   REAL LowerBound() const 
   {
     return std::min(REAL(0.0), this->min()); //*std::min_element(this->begin(), this->end()));
   }

   REAL EvaluatePrimal() const 
   {
     //if(primal_ == no_primal_decision) { primal_ = no_primal_active; }

     if(primal_ < size()) {
       return (*this)[primal_];
     } else if(primal_ == no_primal_active || no_primal_decision) {
       return 0.0;
     } else {
       return std::numeric_limits<REAL>::infinity();
     }
   }

   void init_primal() { primal_ = no_primal_active; }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }
   auto export_variables() { return std::tie( *static_cast<vector<REAL>*>(this) ); }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::vector vars) const
  {
    s.add_at_most_one_constraint(vars.begin(), vars.end());
  }

  /*
  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    sat_literal_vector literals(this->size()+1);
    sat_literal sum_literal;
    load_sat_literals(begin, literals, sum_literal);

    const REAL min_cost = this->min();
    if(min_cost > th) {
      assumptions.push_back(-sum_literal);
    } else {
      if(min_cost < -th) {
        assumptions.push_back(sum_literal);
      }
      for(INDEX i=0; i<this->size(); ++i) {
        if((*this)[i] > min_cost + th) { 
          assumptions.push_back(-literals[i]);
        }
      }
    }
  }
  */

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::vector vars)
  {
    if(s.solution(vars[0])) {
      primal_ = 0; 
    } else if(s.solution(vars[1])) {
      primal_ = 1; 
    } else {
      primal_ = no_primal_active;
    }
  }

private:
   INDEX primal_;
};

// connecting first entry of detection factor with one entry in at_most_one_cell_factor
// left factor is detection factor, right one at_most_one_cell_factor
class at_most_one_cell_message {
public:
  at_most_one_cell_message(const INDEX at_most_one_cell_factor_index) : at_most_one_cell_factor_index_(at_most_one_cell_factor_index) {}

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l.detection += msg;
  }

  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r[at_most_one_cell_factor_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    if(l.detection_active()) {
        r.primal_ = at_most_one_cell_factor_index_;
    }
  }

  // it is more effective to send messages one by one
  template<typename LEFT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToRight_deactivated(const LEFT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    const REAL delta = l.min_detection_cost();
    assert(!std::isnan(delta));

    INDEX no_messages = 0;
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) { ++no_messages; }
    for(; msg_begin!=msg_end; ++msg_begin) {
      (*msg_begin)[0] -= omega/REAL(no_messages)*delta;
    } 
  }

  template<typename LEFT_FACTOR, typename MSG_ITERATOR>
  static REAL send_messages_to_right_improvement(LEFT_FACTOR& l, MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
  {
      const auto diff = l.min_detection_cost();
      const auto left_lb = l.LowerBound();
      const auto new_left_lb = 0.0;

      std::size_t no_messages = 0;
      REAL right_lb = 0.0;
      REAL new_right_lb = 0.0;
      for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) { ++no_messages; }
      for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
          auto* r = (*msg_it)->GetRightFactor()->GetFactor();
          const auto index = (*msg_it).GetMessageOp().at_most_one_cell_factor_index_;
          right_lb += r->LowerBound();
          (*r)[index] += 1.0/REAL(no_messages)*diff;
          new_right_lb += r->LowerBound();
          (*r)[index] -= 1.0/REAL(no_messages)*diff;
      }

      assert(new_left_lb + new_right_lb - (left_lb + right_lb) >= -eps); 
      return new_left_lb + new_right_lb - (left_lb + right_lb);
  }


  template<typename RIGHT_FACTOR, typename MSG>
  void ReceiveRestrictedMessageFromRight(RIGHT_FACTOR& r, MSG& msg) 
  {
    assert(false);
    assert(at_most_one_cell_factor_index_ < r.size());
    //std::cout << "r.primal = " << r.primal_ << " = " << r.no_primal_active << "\n";
    if(r.primal_ == r.no_primal_active || r.primal_ == r.no_primal_decision) { // no element chosen yet
      //make_right_factor_uniform(r, msg);
      const REAL cur_detection_cost = r[at_most_one_cell_factor_index_];
      r[at_most_one_cell_factor_index_] = std::numeric_limits<REAL>::infinity();
      const REAL rest_cost = std::min(REAL(0.0), r.min()); //*std::min_element(r.begin(), r.end()));
      r[at_most_one_cell_factor_index_] = cur_detection_cost;

      msg[0] -= (cur_detection_cost - rest_cost);
    } else if(r.primal_ < r.size()) { // one element already chosen
      if(r.primal_ == at_most_one_cell_factor_index_) {
//        const REAL val_prev = msg.GetLeftFactor()->GetFactor()->detection;
//        msg[0] -= -100000;//std::numeric_limits<REAL>::infinity();
//        const REAL val = msg.GetLeftFactor()->GetFactor()->detection;
//        //assert(msg.GetLeftFactor()->GetFactor()->operator[](0) < -10000);
      } else {
//        //assert(r.primal_ < r.size());
//        const REAL test_val_prev = msg.GetLeftFactor()->GetFactor()->detection;
//        msg[0] -= 100000;//std::numeric_limits<REAL>::infinity();
//        const REAL test_val = msg.GetLeftFactor()->GetFactor()->detection;
//        //assert(test_val > 10000);
//        //assert(msg.GetLeftFactor()->GetFactor()->operator[](0) > 10000);
      }
    } else {
//      assert(r.primal_ == r.primal_infeasible);
//      const REAL test_val_prev = msg.GetLeftFactor()->GetFactor()->detection;
//      msg[0] -= 100000;//std::numeric_limits<REAL>::infinity();
//      const REAL test_val = msg.GetLeftFactor()->GetFactor()->detection;
      //assert(test_val > 10000);
      //assert(msg.GetLeftFactor()->GetFactor()->operator[](0) > 10000); 
    }
  }

  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    //assert(std::distance(msg_begin, msg_end) <= rightFactor.size()); // std::distance does not work currently on the message iterator.
    std::vector<bool> edge_taken(rightFactor.size(), false);
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
        auto& msg = (*msg_it).GetMessageOp();
        const auto idx = msg.at_most_one_cell_factor_index_;
        assert( idx < edge_taken.size() );
        assert(edge_taken[ idx ] == false);
        edge_taken[ idx ] = true;
    }
    assert(std::count(edge_taken.begin(), edge_taken.end(), false) == 0);

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = 0.0;
    
    const auto right_factor_size = rightFactor.size();
    for(INDEX i=0; i<right_factor_size; ++i) {
      const REAL val = rightFactor[i];
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    REAL set_to_cost = std::min(smallest_not_taken, second_smallest_taken);

    for(; msg_begin!=msg_end; ++msg_begin) {
      (*msg_begin)[0] -= omega*( rightFactor[ (*msg_begin).GetMessageOp().at_most_one_cell_factor_index_ ] - set_to_cost );
    }
  }

  template<typename RIGHT_FACTOR, typename MSG_ITERATOR>
  static REAL send_messages_to_left_improvement(RIGHT_FACTOR& r, MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
  {
      assert(false);
      return 0.0; 
  }


  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
  {
    // make cost of detection same as cost of non-detection
    msg[0] -= omega*l.min_detection_cost();
  }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    REAL send_message_to_right_improvement(LEFT_FACTOR& l, RIGHT_FACTOR& r) const
    {
        const auto left_lb = l.LowerBound();
        const auto new_left_lb = 0.0;

        const auto diff = l.min_detection_cost();

        const auto right_lb = r.LowerBound();
        const auto prev_val = r[at_most_one_cell_factor_index_];
        r[at_most_one_cell_factor_index_] += diff;
        const auto new_right_lb = r.LowerBound();
        r[at_most_one_cell_factor_index_] = prev_val;

        assert(new_left_lb + new_right_lb - (left_lb + right_lb) >= -eps);
        return new_left_lb + new_right_lb - (left_lb + right_lb);
    }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
  {
    const REAL cur_detection_cost = r[at_most_one_cell_factor_index_];
    r[at_most_one_cell_factor_index_] = std::numeric_limits<REAL>::infinity();
    const REAL rest_cost = std::min(REAL(0.0), r.min());
    r[at_most_one_cell_factor_index_] = cur_detection_cost;

    msg[0] -= omega*(cur_detection_cost - rest_cost);
  }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    REAL send_message_to_left_improvement(LEFT_FACTOR& l, RIGHT_FACTOR& r) const
    {
        const auto cur_detection_cost = r[at_most_one_cell_factor_index_];
        r[at_most_one_cell_factor_index_] = std::numeric_limits<REAL>::infinity();
        const auto rest_cost = r.LowerBound();
        r[at_most_one_cell_factor_index_] = cur_detection_cost;

        const auto diff_right = rest_cost;
        const auto left_lb = l.LowerBound();
        const auto new_left_lb = std::min(0.0, l.min_incoming() + r.detection + (cur_detection_cost - rest_cost) + r.min_outgoing());

        return diff_right + new_left_lb - left_lb;
    }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    //assert(false); 
    if(r.primal_ == r.primal_infeasible) {
      //std::cout << "exclusion constraints not satisfied!\n";
      return false;
    }
    if(l.detection_active()) {
      return (r.primal_ == at_most_one_cell_factor_index_);
    } else {
      return (r.primal_ != at_most_one_cell_factor_index_);
    }
  }

  // without division distance
  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable detection_var, typename SOLVER::vector incoming_vars, typename SOLVER::vector outgoing_vars,
      RIGHT_FACTOR& r, typename SOLVER::vector at_most_one_detection_vars) const
  {
    s.make_equal(detection_var, at_most_one_detection_vars[at_most_one_cell_factor_index_]);
  }
  // with division distance
  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::vector detection_vars, typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars, typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars,
      RIGHT_FACTOR& r, typename SOLVER::vector at_most_one_detection_vars) const
  {
    auto detection_active = s.add_at_most_one_constraint(detection_vars.begin(), detection_vars.end());
    s.make_equal(detection_active, at_most_one_detection_vars[at_most_one_cell_factor_index_]);
  }
  
  // for flow conservation factor
  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable detection_var, typename SOLVER::vector left_edges, 
      RIGHT_FACTOR& r, typename SOLVER::vector at_most_one_detection_vars) const
  {
    s.make_equal(detection_var, at_most_one_detection_vars[at_most_one_cell_factor_index_]);
  }
  // with division distance

private:
   const INDEX at_most_one_cell_factor_index_;
};

// first entry: lower exit transition. second entry: upper transition messages (excluding exit transition)
// x_1 = 1 implies x_2 = 0
// equivalent to x_1 + x_2 <= 1
class exit_constraint_factor : public array<REAL,2> {
  template<exit_constraint_position POSITION> friend class exit_constraint_message;
public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_primal_active = std::numeric_limits<INDEX>::max()-1;
  constexpr static INDEX primal_infeasible = std::numeric_limits<INDEX>::max()-2; 

  using repam_type = array<REAL,2>;
  REAL LowerBound() const 
  {
    return std::min(REAL(0.0), std::min( (*this)[0], (*this)[1] )); 
  }
  REAL EvaluatePrimal() const 
  {
    if(primal_ < this->size()) {
      return (*this)[primal_];
    } else if(primal_ == no_primal_active) {
      return 0.0;
    } else {
      return std::numeric_limits<REAL>::infinity();
    }
  }

  void init_primal() { primal_ = no_primal_decision; }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( binary_data<REAL>(&(*this)[0], 2) ); }
  auto export_variables() { return std::tie( *static_cast<repam_type*>(this) ); }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::vector vars) const
  {
    s.add_at_most_one_constraint(vars.begin(), vars.end());
  }

  /*
  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    sat_literal_vector literals(this->size()+1);
    load_sat_literals(begin, literals);

    // do zrobienia: if reparametrization is < -th, then disallow no edge to be taken!
    for(INDEX i=0; i<this->size(); ++i) {
      if((*this)[i] > th) { 
        assumptions.push_back(-literals[i]); 
      }
    }
  }
  */

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::vector vars)
  {
    if(s.no_active(vars.begin(), vars.end()) > 0) {
      primal_ = s.first_active(vars);
    } else {
      primal_ = no_primal_active;
    }
  }

private:
  INDEX primal_;
};

// left is detection_factor, right is exit_constraint_factor
template<exit_constraint_position POSITION>
class exit_constraint_message {
public:
  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    if(POSITION == exit_constraint_position::lower) {
      const INDEX i = l.outgoing.size()-1;
      l.update_outgoing(i, msg);
    } else {
      const auto outgoing_size = l.outgoing.size();
      for(INDEX i=0; i<outgoing_size-1; ++i) {
        l.update_outgoing(i, msg);
      }
    }
  }

  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    if(POSITION == exit_constraint_position::lower) {
      r[0] += msg;
    } else {
      r[1] += msg;
    }
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    //assert(false);
    if(POSITION == exit_constraint_position::lower) {
      if(l.outgoing_edge_ == l.outgoing.size() - 1) {
        r.primal_ = 0;
      } 
    } else {

    }
  }

  template<typename LEFT_FACTOR, typename MSG>
  void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
  { 
    REAL detection_incoming_cost = l.detection;
    detection_incoming_cost += l.incoming.min(); //*std::min_element(l.incoming.begin(), l.incoming.end());
    REAL outgoing_cost = std::numeric_limits<REAL>::infinity();
    assert(l.outgoing.size() > 0);
    if(l.outgoing.size() > 1) {
      outgoing_cost = l.outgoing.min(); //*std::min_element(l.outgoing.begin(), l.outgoing.end()-1);
    } 

    if(POSITION == exit_constraint_position::lower) {
      msg[0] -= omega*( l.outgoing[ l.outgoing.size()-1 ] - std::min( -detection_incoming_cost, outgoing_cost ) );
    } else {
      msg[0] -= omega*( outgoing_cost - std::min( -detection_incoming_cost, l.outgoing[ l.outgoing.size()-1 ] ) );
    }
  }


  template<typename LEFT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToRight(const LEFT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    REAL detection_incoming_cost = l.detection;
    detection_incoming_cost += l.incoming.min(); //*std::min_element(l.incoming.begin(), l.incoming.end());
    REAL outgoing_cost = std::numeric_limits<REAL>::infinity();
    assert(l.outgoing.size() > 0);
    if(l.outgoing.size() > 1) {
      outgoing_cost = l.outgoing.min(); //*std::min_element(l.outgoing.begin(), l.outgoing.end()-1);
    } 

    INDEX no_messages = 0;
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) { no_messages++; }

    if(POSITION == exit_constraint_position::lower) {
      const REAL delta = l.outgoing[ l.outgoing.size()-1 ] - std::min( -detection_incoming_cost, outgoing_cost );
      assert(!std::isnan(delta));
      for(; msg_begin!=msg_end; ++msg_begin) { (*msg_begin)[0] -=  omega/REAL(no_messages)*delta; }
    } else {
      const REAL delta = outgoing_cost - std::min( -detection_incoming_cost, l.outgoing[ l.outgoing.size()-1 ] );
      assert(!std::isnan(delta));
      for(; msg_begin!=msg_end; ++msg_begin) { (*msg_begin)[0] -=  omega/REAL(no_messages)*delta; }
    }
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega) 
  {
    assert(r.size() == 2);
    if(POSITION == exit_constraint_position::lower) {
      msg[0] -= omega*(r[0] - std::min(0.0, r[1]));
    } else {
      msg[0] -= omega*(r[1] - std::min(0.0, r[0]));
    }
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    if(POSITION == exit_constraint_position::lower) {
      if(r.primal_ == 0) {
        return l.outgoing_edge_ == l.outgoing.size()-1;
      } else {
        return l.outgoing_edge_ < l.outgoing.size()-1 || l.outgoing_edge_ == l.no_edge_taken;
      }
    } else {
      if(r.primal_ == 1) {
        return l.outgoing_edge_ < l.outgoing.size()-1;
      } else {
        return l.outgoing_edge_ == l.outgoing.size()-1 || l.outgoing_edge_ == l.no_edge_taken;
      }
    }
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable detection_var, typename SOLVER::vector incoming_vars, typename SOLVER::vector outgoing_vars,
      RIGHT_FACTOR& r, typename SOLVER::vector exit_vars) const
  {
    if(POSITION == exit_constraint_position::lower) {
      s.make_equal(outgoing_vars[outgoing_vars.size()-1], exit_vars[0]);
    } else {
      // get indicator variable for sum of outgoing edges without last one
      auto c = s.add_at_most_one_constraint(outgoing_vars.begin(), outgoing_vars.end()-1);
      s.make_equal(c, exit_vars[1]);
    }
  }
};


// detection factor disallowing division when it occured less than division_distance timesteps before.
// degenerates to ordinary detection_factor when division_distance = 1.
// possible rename _dd to _with_division_distance
class detection_factor_dd {
  friend class at_most_one_cell_message;

public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_edge_taken = std::numeric_limits<INDEX>::max()-1;

  detection_factor_dd(
      const INDEX _no_incoming_transition_edges, const INDEX _no_incoming_division_edges, const INDEX _no_outgoing_transition_edges, const INDEX _no_outgoing_division_edges, 
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, const INDEX _division_distance
      )
    :
      detection(_division_distance, detection_cost), // to do: detection cost need not be a vector, scalar is enough!
      incoming_division(_no_incoming_division_edges+1, std::numeric_limits<REAL>::infinity()),
      incoming_transition(_no_incoming_transition_edges+1, _division_distance-1, std::numeric_limits<REAL>::infinity()),
      outgoing_transition(_no_outgoing_transition_edges+1, _division_distance, std::numeric_limits<REAL>::infinity()),
      outgoing_division(_no_outgoing_division_edges+1, std::numeric_limits<REAL>::infinity())
  {
    for(INDEX i=0; i<detection.size(); ++i) { 
      assert(detection[i] == detection_cost); 
    }
    assert(_division_distance >= 2); 
    incoming_division.back() = appearance_cost;
    outgoing_division.back() = disappearance_cost;
    for(INDEX i=0; i<division_distance()-1; ++i) {
      incoming_transition(no_incoming_transition_edges()-1, i) = appearance_cost;
    }
    for(INDEX i=0; i<division_distance(); ++i) {
      outgoing_transition(no_outgoing_transition_edges()-1, i) = disappearance_cost;
    }
  }

  INDEX division_distance() const { return detection.size(); }
  INDEX no_incoming_division_edges() const { return incoming_division.size(); }
  INDEX no_incoming_transition_edges() const { return incoming_transition.dim1(); }
  INDEX no_outgoing_transition_edges() const { return outgoing_transition.dim1(); }
  INDEX no_outgoing_division_edges() const { return outgoing_division.size(); }

  INDEX size() const
  {
    //assert(false); 
    // note: appearance and disappearance are stored separately
    return 0;
    //return division_distance() + no_incoming_transition_edges()*(division_distance()-1) + no_incoming_division_edges() + no_outgoing_transition_edges()*division_edges() + no_outgoing_division_edges();
  }

  void set_incoming_transition_cost(const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX edge_index, const REAL cost)
  {
    assert(edge_index < this->no_incoming_transition_edges()-1);
    for(INDEX i=0; i<division_distance()-1; ++i) {
      assert(incoming_transition(edge_index, i) == std::numeric_limits<REAL>::infinity());
      incoming_transition(edge_index, i) = cost;
    }
  }

  void set_outgoing_transition_cost(const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, const INDEX edge_index, const REAL cost)
  {
    assert(edge_index < this->no_outgoing_transition_edges()-1);
    for(INDEX i=0; i<division_distance(); ++i) {
      assert(outgoing_transition(edge_index, i) == std::numeric_limits<REAL>::infinity());
      outgoing_transition(edge_index, i) = cost;
    }
  }

  void set_incoming_division_cost(const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX edge_index, const REAL cost)
  {
    assert(edge_index < this->no_incoming_division_edges()-1);
    assert(incoming_division[edge_index] == std::numeric_limits<REAL>::infinity());
    incoming_division[edge_index] = cost;
  }

  void set_outgoing_division_cost(const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, const INDEX edge_index, const REAL cost)
  {
    assert(edge_index < this->no_outgoing_division_edges()-1);
    assert(outgoing_division[edge_index] == std::numeric_limits<REAL>::infinity());
    outgoing_division[edge_index] = cost;
  }

  REAL min_detection_cost() const
  {
    auto outgoing_transition_min = outgoing_transition.min2();
    assert(outgoing_transition_min.size() == division_distance());
    auto incoming_transition_min = incoming_transition.min2();
    assert(incoming_transition_min.size() == division_distance()-1);

    const REAL c_first = detection[0] + incoming_division.min() + outgoing_transition_min[0];
    REAL cost = c_first;

    for(INDEX i=1; i<division_distance()-1; ++i) {
      const REAL c = detection[i] + incoming_transition_min[i-1] + outgoing_transition_min[i];
      cost = std::min(cost, c);
    }

    const INDEX last = division_distance()-1;
    const REAL c_last = detection[last] + incoming_transition_min[last-1] + std::min(outgoing_transition_min[last], outgoing_division.min());

    cost = std::min(cost, c_last);
    return cost;
  }

  REAL LowerBound() const
  {
    // assert that all costs have been set ///////////////////////////////////////////////////////////////////////////////////
    for(INDEX i=0; i<detection.size(); ++i) { assert(detection[i] < std::numeric_limits<REAL>::infinity()); }
    for(INDEX i=0; i<incoming_division.size(); ++i) { assert(incoming_division[i] < std::numeric_limits<REAL>::infinity()); }
    for(INDEX i=0; i<outgoing_division.size(); ++i) { assert(outgoing_division[i] < std::numeric_limits<REAL>::infinity()); }
    for(INDEX i=0; i<incoming_transition.dim1(); ++i) {
      for(INDEX j=0; j<incoming_transition.dim2(); ++j) {
        assert(incoming_transition(i,j) < std::numeric_limits<REAL>::infinity()); 
      }
    }
    for(INDEX i=0; i<outgoing_transition.dim1(); ++i) {
      for(INDEX j=0; j<outgoing_transition.dim2(); ++j) {
        assert(outgoing_transition(i,j) < std::numeric_limits<REAL>::infinity()); 
      }
    }


    //std::cout << detection.min() << ", " << incoming_transition.min() << ", " << outgoing_transition.min() << ", " << incoming_division.min() << ", " << outgoing_division.min() << "\n";
    //std::cout << min_detection_cost() << "\n";
    //exit(1);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return std::min(min_detection_cost(), REAL(0.0)); 
  }

  REAL EvaluatePrimal() const
  {
    if(primal().division == no_edge_taken) {
      return 0.0;
    } else if(primal().division < division_distance()) {
      if(primal().division == 0) {
        return detection[0] + incoming_division[primal().incoming_edge] + outgoing_transition(primal().outgoing_edge, 0);
      } else if(primal().division < division_distance()-1) {
        return detection[primal().division] + incoming_transition(primal().incoming_edge, primal().division-1) + outgoing_transition(primal().outgoing_edge, primal().division);
      } else {
        assert(primal().division == division_distance()-1);
        const REAL incoming_cost = incoming_transition(primal().incoming_edge, primal().division-1);
        const REAL outgoing_cost = primal().outgoing_division ? outgoing_division[primal().outgoing_edge] : outgoing_transition(primal().outgoing_edge, primal().division);
        return detection[primal().division] + incoming_cost + outgoing_cost;
      }
    } else {
      return std::numeric_limits<REAL>::infinity();
    }
  }

  //void init_primal() { primal().incoming_edge = no_primal_decision; primal().outgoing_edge = no_primal_decision; primal().division = no_primal_decision; }
  void init_primal() { primal().incoming_edge = no_edge_taken; primal().outgoing_edge = no_edge_taken; primal().division = no_edge_taken; primal().outgoing_division = false; }
  bool detection_active() const { return primal().division < division_distance(); }

  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal().incoming_edge, primal().outgoing_edge, primal().division, primal().outgoing_division ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( detection, incoming_division, incoming_transition, outgoing_transition, outgoing_division ); }
  auto export_variables() { return std::tie( detection, incoming_division, incoming_transition, outgoing_transition, outgoing_division ); }

  // reparametrization
  vector<REAL> detection;
  vector<REAL> incoming_division;
  matrix<REAL> incoming_transition; // (edge no, division distance)
  matrix<REAL> outgoing_transition; // (edge no, division distance)
  vector<REAL> outgoing_division;

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::vector detection_vars, 
      typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars,
      typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars ) const
  {
    auto detection_active = s.add_at_most_one_constraint(detection_vars.begin(), detection_vars.end());

    // detection var must be equal to incoming and outgoing var
    s.make_equal(detection_vars[0], s.add_at_most_one_constraint(incoming_division_vars.begin(), incoming_division_vars.end()));
    s.make_equal(detection_vars[0], s.add_at_most_one_constraint(outgoing_transition_vars.slice2(0).begin(), outgoing_transition_vars.slice2(0).end()));

    for(INDEX i=1; i<detection.size()-1; ++i) {
      auto incoming_slice = incoming_transition_vars.slice2(i-1);
      auto outgoing_slice = outgoing_transition_vars.slice2(i);
      s.make_equal( detection_vars[i], s.add_at_most_one_constraint(incoming_slice.begin(), incoming_slice.end()));
      s.make_equal( detection_vars[i], s.add_at_most_one_constraint(outgoing_slice.begin(), outgoing_slice.end()));
    }

    auto last_incoming_slice = incoming_transition_vars.slice2( incoming_transition_vars.dim2()-1 );
    s.make_equal( detection_vars.back(), s.add_at_most_one_constraint(last_incoming_slice.begin(), last_incoming_slice.end()) );

    auto last_outgoing_transition_slice = outgoing_transition_vars.slice2( outgoing_transition_vars.dim2()-1 );
    auto last_outgoing_transition_sum = s.add_at_most_one_constraint( last_outgoing_transition_slice.begin(), last_outgoing_transition_slice.end() );
    auto outgoing_division_literal = s.add_at_most_one_constraint( outgoing_division_vars.begin(), outgoing_division_vars.end() );
    std::array<typename SOLVER::variable,2> last_outgoing_vars({last_outgoing_transition_sum, outgoing_division_literal}); 
    auto last_outgoing_var = s.add_at_most_one_constraint(last_outgoing_vars.begin(), last_outgoing_vars.end());
  }

  /*
  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var first_detection_var) const
  {
    sat_literal detection_active_literal;
    sat_literal_vector detection_literals(detection);
    sat_literal_vector incoming_division_literals(incoming_division);
    sat_literal_matrix incoming_transition_literals(incoming_transition);
    sat_literal_matrix outgoing_transition_literals(outgoing_transition);
    sat_literal_vector outgoing_division_literals(outgoing_division);
    load_sat_literals(detection_active_literal, detection_literals, incoming_division_literals, incoming_transition_literals, outgoing_transition_literals, outgoing_division_literals); 

    auto outgoing_transition_min = outgoing_transition.min2();
    assert(outgoing_transition_min.size() == division_distance());
    auto incoming_transition_min = incoming_transition.min2();
    assert(incoming_transition_min.size() == division_distance()-1); 

    const REAL min_det_cost = min_detection_cost();
    //std::cout << "\nreduce sat: " << first_detection_var << ";" << min_det_cost <<  "," << th << "\n";
    //std::cout << "no incoming edges: " << no_incoming_transition_edges() << "\n";
    //std::cout << "no outgoing edges: " << no_outgoing_transition_edges() << "\n";

//    if(no_incoming_transition_edges() == 2) {
//      std::cout << "fix first incoming transition in second timestep\n";
//      assumptions.push_back(to_literal(first_incoming_transition));
//    }
//    if(no_outgoing_transition_edges() == 2) {
//      std::cout << "fix first outgoing transition in first timestep\n";
//      assumptions.push_back(to_literal(first_outgoing_transition));
//    }
//    return;

    for(INDEX i=0; i<no_incoming_division_edges(); ++i) {
      //std::cout << incoming_division[i] << ",";
    }
    //std::cout << "\n\n";
    for(INDEX t=0; t<division_distance()-1; ++t) {
      for(INDEX i=0; i<no_incoming_transition_edges(); ++i) {
        //std::cout << incoming_transition(i,t) << ",";
      }
      //std::cout << "\n";
    }
    //std::cout << "\n";
    for(INDEX t=0; t<division_distance(); ++t) {
      for(INDEX i=0; i<no_outgoing_transition_edges(); ++i) {
        //std::cout << outgoing_transition(i,t) << ",";
      }
      //std::cout << "\n";
    }
    //std::cout << "\n";
    for(INDEX i=0; i<no_outgoing_division_edges(); ++i) {
      //std::cout << outgoing_division[i] << ",";
    }
    //std::cout << "\n\n";

    if(min_det_cost > th) {
      assumptions.push_back(-detection_active_literal);
      //for(INDEX j=0; j<incoming_division.size(); ++j) {
      //  assumptions.push_back(-to_literal(begin + 1 + j));
      //} 
    } else {
      if(min_det_cost <= -th) {
        //std::cout << "detection must be on\n";
        assumptions.push_back(detection_active_literal);
      }
      const REAL incoming_division_min = incoming_division.min();
      const REAL det_0_cost = detection[0] + incoming_division_min + outgoing_transition_min[0];
      if(det_0_cost <= min_det_cost + th) {
        //std::cout << "first detection on\n";
        // go over incoming division edges
        for(INDEX i=0; i<incoming_division.size(); ++i) {
          if(incoming_division[i] > incoming_division_min + th) {
            assumptions.push_back(-incoming_division_literals[i]);
          }
        }
        for(INDEX i=0; i<no_outgoing_transition_edges(); ++i) {
          if(outgoing_transition(i,0) > outgoing_transition_min[0] + th) {
            assumptions.push_back(-outgoing_transition_literals(i,0));
            //std::cout << "forbid outgoing transition " << i << ",0" << "\n";
          } 
        }
      } else {
        assumptions.push_back(-detection_literals[0]);
      }

      for(INDEX t=1; t<division_distance()-1; ++t) {
        if(detection[t] + incoming_transition_min[t-1] + outgoing_transition_min[t] <= min_det_cost + th) {
          for(INDEX i=0; i<no_incoming_transition_edges(); ++i) {
            if(incoming_transition(i,t-1) > incoming_transition_min[t-1] + th) {
              assumptions.push_back(-incoming_transition_literals(i,t-1));
            }
          }
          for(INDEX i=0; i<no_outgoing_transition_edges(); ++i) {
            if(outgoing_transition(i,t) > outgoing_transition_min[t] + th) {
              assumptions.push_back(-outgoing_transition_literals(i,t));
            }
          } 
        } else {
          assumptions.push_back(-detection_literals[t]);
        } 
      }

      const INDEX last = division_distance()-1;
      if(detection[last] + incoming_transition_min[last-1] + std::min(outgoing_transition_min[last], outgoing_division.min()) <= min_det_cost + th) {
//        std::cout << "last detection on\n";
        for(INDEX i=0; i<no_incoming_transition_edges(); ++i) {
          if(incoming_transition(i,last-1) > incoming_transition_min[last-1] + th) {
            assumptions.push_back(-incoming_transition_literals(i,last-1));
          }
        }
        const REAL outgoing_min = std::min(outgoing_transition_min[last], outgoing_division.min());
        for(INDEX i=0; i<no_outgoing_transition_edges(); ++i) {
          if(outgoing_transition(i,last) > outgoing_min + th) {
            assumptions.push_back(-outgoing_transition_literals(i,last));
          }
        } 
        for(INDEX i=0; i<outgoing_division.size(); ++i) {
          if(outgoing_division[i] > outgoing_min + th) {
            assumptions.push_back(-outgoing_division_literals[i]);
          }
        }
      } else {
        assumptions.push_back(-detection_literals[last]);
      }
    }
    //std::cout << "\n";
  }
*/

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::vector detection_vars, 
      typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars,
      typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars )
  {
    primal().division = no_edge_taken; 
    primal().incoming_edge = no_edge_taken;
    primal().outgoing_edge = no_edge_taken; 
    primal().outgoing_division = false;

    if(!s.no_active(detection_vars.begin(), detection_vars.end())) { 
      //std::cout << "no detection\n";
      return; 
    }

    assert(false); // check whether slices are done correctly below.
    const INDEX timestep = s.first_active(detection_vars);
    primal().division = timestep;
    if(primal().division == 0) {
      assert(primal().incoming_edge == no_edge_taken);
      primal().incoming_edge = s.first_active(incoming_division_vars);
      auto slice = outgoing_transition_vars.slice_left(0);
      primal().outgoing_edge = s.first_active(slice.begin(), slice.end());
    } else {
      auto incoming_slice = incoming_transition_vars.slice_left(timestep);
      primal().incoming_edge = s.first_active(incoming_slice.begin(), incoming_slice.end());
      auto outgoing_slice = outgoing_transition_vars.slice_left(timestep);
      if(primal().division < division_distance()-1) {
        primal().outgoing_edge = s.first_active(outgoing_slice.begin(), outgoing_slice.end());
      } else {
        if(s.no_active(outgoing_division_vars.begin(), outgoing_division_vars.end()) > 0) {
          primal().outgoing_edge = s.first_active(outgoing_division_vars);
          primal().outgoing_division = true;
        } else {
          primal().outgoing_edge = s.first_active(outgoing_slice.begin(), outgoing_slice.end()); 
        }
      }
    }
  }
  bool primal_valid() const
  {
    if(primal().division >= division_distance()) { return false; }
    if(primal().outgoing_division) {
      if(primal().outgoing_edge >= no_outgoing_division_edges()) { return false; }
    } else {
      if(primal().outgoing_edge >= no_outgoing_transition_edges()) { return false; }
    }
    if(primal().division == 0) {
      if(primal().incoming_edge >= no_incoming_division_edges()) { return false; }
    } else { 
      if(primal().incoming_edge >= no_incoming_transition_edges()) { return false; }
    }
    //std::cout << "\n";
    //std::cout << primal().division << ";" << primal().incoming_edge << ";" << primal().outgoing_division << "," << primal().outgoing_edge << "\n";
    return true;
  }

  INDEX division() const { return primal().division; }

  struct primal_struct {
    INDEX incoming_edge, outgoing_edge, division;
    bool outgoing_division;
  };

  const primal_struct& primal() const { return p_; }
  primal_struct& primal() { return p_; }

  REAL min_incoming_division_marginal_diff(const INDEX incoming_edge_index)
  {
    assert(incoming_edge_index < no_incoming_division_edges());

    auto outgoing_transition_min = outgoing_transition.min2();
    assert(outgoing_transition_min.size() == division_distance());
    auto incoming_transition_min = incoming_transition.min2();
    assert(incoming_transition_min.size() == division_distance()-1);

    const REAL cost_1 = detection[0] + outgoing_transition_min[0] + incoming_division[ incoming_edge_index ];

    const REAL incoming_val = incoming_division[incoming_edge_index];
    incoming_division[incoming_edge_index] = std::numeric_limits<REAL>::infinity();
    const REAL min_incoming_val = incoming_division.min();
    incoming_division[incoming_edge_index] = incoming_val;

    REAL cost_0 = detection[0] + outgoing_transition_min[0] + min_incoming_val;

    for(INDEX i=1; i<division_distance()-1; ++i) {
      const REAL c = detection[i] + incoming_transition_min[i-1] + outgoing_transition_min[i];
      cost_0 = std::min(c, cost_0);
    }

    const INDEX last = division_distance()-1;
    const REAL c_last = detection[last] + incoming_transition_min[last-1] + std::min(outgoing_transition_min[last], outgoing_division.min());
    cost_0 = std::min(c_last, cost_0);

    cost_0 = std::min(REAL(0.0), cost_0); 

    return cost_1 - cost_0; 
  }

  vector<REAL> min_incoming_transition_marginal_diff(const INDEX incoming_edge_index)
  {
    assert(incoming_edge_index < no_incoming_transition_edges());

    vector<REAL> original_edge_val(incoming_transition.dim2());
    assert(original_edge_val.size() == division_distance()-1);
    for(INDEX i=0; i<original_edge_val.size(); ++i) {
      original_edge_val[i] = incoming_transition(incoming_edge_index, i);
      incoming_transition(incoming_edge_index, i) = std::numeric_limits<REAL>::infinity();
    }
    vector<REAL> incoming_transition_min = incoming_transition.min2();
    assert(incoming_transition_min.size() == division_distance()-1);
    for(INDEX i=0; i<original_edge_val.size(); ++i) {
      incoming_transition(incoming_edge_index, i) = original_edge_val[i];
    }

    vector<REAL> outgoing_transition_min = outgoing_transition.min2();
    assert(outgoing_transition_min.size() == division_distance());
    const REAL outgoing_division_min = outgoing_division.min();

    const INDEX last = division_distance()-1;

    REAL cost_not_taken = 0.0;
    std::array<REAL,2> smallest_taken = {std::numeric_limits<REAL>::infinity(),std::numeric_limits<REAL>::infinity()};
    auto update_taken = [&smallest_taken](const REAL x) {
      const REAL min = std::min(smallest_taken[0], x);
      const REAL max = std::max(smallest_taken[0], x);
      smallest_taken[0] = min;
      smallest_taken[1] = std::min(max, smallest_taken[1]); 
    };
    auto update_not_taken = [&cost_not_taken](const REAL x) { cost_not_taken = std::min(cost_not_taken, x); }; 

    update_not_taken( detection[0] + incoming_division.min() + outgoing_transition_min[0] );
    for(INDEX i=1; i<last; ++i) {
      update_not_taken( detection[i] + incoming_transition_min[i-1] + outgoing_transition_min[i] );
    }
    update_not_taken( detection[last] + incoming_transition_min[last-1] + std::min(outgoing_transition_min[last], outgoing_division_min) );

    vector<REAL> cost_diff(division_distance()-1); 
    for(INDEX i=1; i<last; ++i) {
      const REAL val = detection[i] + original_edge_val[i-1] + outgoing_transition_min[i];
      update_taken(val);
      cost_diff[i-1] = val;
    }
    {
      const REAL last_val = detection[last] + original_edge_val[last-1] + std::min(outgoing_transition_min[last], outgoing_division_min);
      update_taken(last_val);
      cost_diff[last-1] = last_val;
    }

    const REAL set_to_cost = std::min(cost_not_taken, smallest_taken[1]);
    //const REAL set_to_cost = std::min(cost_not_taken, smallest_taken[0]);
    for(INDEX i=0; i<cost_diff.size(); ++i) {
      cost_diff[i] -= set_to_cost;
    }

    return std::move(cost_diff); 
  }

  vector<REAL> min_outgoing_transition_marginal_diff(const INDEX outgoing_edge_index)
  {
    assert(outgoing_edge_index < no_outgoing_transition_edges());

    vector<REAL> original_edge_val(outgoing_transition.dim2());
    assert(original_edge_val.size() == division_distance());
    for(INDEX i=0; i<division_distance(); ++i) {
      original_edge_val[i] = outgoing_transition(outgoing_edge_index, i);
      outgoing_transition(outgoing_edge_index, i) = std::numeric_limits<REAL>::infinity();
    }
    vector<REAL> outgoing_transition_min = outgoing_transition.min2();
    assert(outgoing_transition_min.size() == division_distance());
    for(INDEX i=0; i<division_distance(); ++i) {
      outgoing_transition(outgoing_edge_index, i) = original_edge_val[i];
    }

    auto incoming_transition_min = incoming_transition.min2();
    assert(incoming_transition_min.size() == division_distance()-1);

    const REAL incoming_division_min = incoming_division.min();

    const INDEX last = division_distance()-1;

    // minimum over all assignments not taken
    REAL cost_not_taken = 0.0;
    std::array<REAL,2> smallest_taken = {std::numeric_limits<REAL>::infinity(),std::numeric_limits<REAL>::infinity()};
    auto update_taken = [&smallest_taken](const REAL x) {
      const REAL min = std::min(smallest_taken[0], x);
      const REAL max = std::max(smallest_taken[0], x);
      smallest_taken[0] = min;
      smallest_taken[1] = std::min(max, smallest_taken[1]); 
    };
    auto update_not_taken = [&cost_not_taken](const REAL x) { cost_not_taken = std::min(cost_not_taken, x); }; 

    update_not_taken( detection[0] + incoming_division_min + outgoing_transition_min[0] );
    for(INDEX i=1; i<division_distance()-1; ++i) {
      update_not_taken( detection[i] + incoming_transition_min[i-1] + outgoing_transition_min[i] );
    }
    update_not_taken( detection[last] + incoming_transition_min[last-1] + std::min(outgoing_transition_min[last], outgoing_division.min()) );

    vector<REAL> cost_diff(division_distance()-1);

    cost_diff[0] = detection[0] + incoming_division_min + original_edge_val[0];
    for(INDEX i=1; i<last-1; ++i) {
      const REAL val = detection[i] + incoming_transition_min[i-1] + original_edge_val[i];
      update_taken(val);
      cost_diff[i] = val;
    }
    {
      const REAL second_last_val = detection[last-1] + incoming_transition_min[last-2] + original_edge_val[last-1];
      const REAL last_val =  detection[last] + incoming_transition_min[last-1] + original_edge_val[last];
      const REAL val = std::min(second_last_val, last_val);
      update_taken(val);
      cost_diff[last-1] = val;
    }

    const REAL set_to_cost = std::min(cost_not_taken, smallest_taken[1]);
    //const REAL set_to_cost = std::min(cost_not_taken, smallest_taken[1]);
    for(INDEX i=0; i<cost_diff.size(); ++i) {
      cost_diff[i] -= set_to_cost;
    }

    return std::move(cost_diff);
  }

  REAL min_outgoing_division_marginal_diff(const INDEX outgoing_edge_index)
  {
    assert(outgoing_edge_index < no_outgoing_division_edges());

    auto outgoing_transition_min = outgoing_transition.min2();
    assert(outgoing_transition_min.size() == division_distance());
    auto incoming_transition_min = incoming_transition.min2();
    assert(incoming_transition_min.size() == division_distance()-1);

    const REAL outgoing_val = outgoing_division[outgoing_edge_index];
    outgoing_division[outgoing_edge_index] = std::numeric_limits<REAL>::infinity();
    const REAL min_outgoing_val = outgoing_division.min();
    outgoing_division[outgoing_edge_index] = outgoing_val; 

    REAL cost_0 = detection[0] + incoming_division.min() + outgoing_transition_min[0];

    for(INDEX i=1; i<division_distance()-1; ++i) {
      const REAL c = detection[i] + incoming_transition_min[i-1] + outgoing_transition_min[i];
      cost_0 = std::min(c, cost_0);
    }

    const INDEX last = division_distance()-1;
    const REAL c_last = detection[last] + incoming_transition_min[last-1] + std::min(outgoing_transition_min[last], min_outgoing_val);
    cost_0 = std::min(c_last, cost_0);

    cost_0 = std::min(REAL(0.0), cost_0); 

    const REAL cost_1 = detection[last] + incoming_transition_min[last-1] + outgoing_val;

    return cost_1 - cost_0; 
  }

private:

  // primal
  primal_struct p_;

  //INDEX incoming_edge_, outgoing_edge_, division_;
  bool outgoing_division_;
};

// is simplex constraint with <= instead of ==. Put at_most_one_cell_factor and this one in simplex_factor.hxx
class mapping_edge_factor_dd : public vector<REAL> {
public:
  static constexpr INDEX no_edge_taken = std::numeric_limits<INDEX>::max()-1;

  //mapping_edge_factor_dd(const INDEX _size) : vector<REAL>(_size,0.0) {}
  mapping_edge_factor_dd(const REAL cost, const INDEX _size) : vector<REAL>(_size-1,cost) { assert(_size > 1); }

  REAL LowerBound() const
  { 
    return std::min(REAL(0.0), this->min());
  }

  REAL EvaluatePrimal() const
  {
    if(primal_ < size()) {
      return (*this)[primal_];
    } else if(primal_ == no_edge_taken) {
      return 0.0;
    } else {
      return std::numeric_limits<REAL>::infinity();
    }
  }
       
  void init_primal() { primal_ = no_edge_taken; }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }
  auto export_variables() { return std::tie( *static_cast<vector<REAL>*>(this) ); }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::vector vars) const
  {
    assert(false);
  }

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::vector vars)
  {
    assert(false);
  }

private: 
  INDEX primal_;
};

class division_edge_factor_dd {
public:
  division_edge_factor_dd(const INDEX _cost) : cost(_cost) {}
  REAL cost;
  REAL LowerBound() const { return std::min(cost, REAL(0.0)); }
  REAL EvaluatePrimal() const {
    if(primal_) {
      return cost;
    } else {
      return 0.0;
    }
  }
  INDEX size() const { return 1; }

  void init_primal() {}
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cost ); }
  auto export_variables() { return std::tie( cost ); }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::variable var) const
  {
    assert(false);
  }

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::variable var)
  {
    assert(false);
  }

private:
  bool primal_;

};

// left is mapping_edge_factor_dd, right is detection_factor_dd
class cell_incoming_mapping_edge_detection_factor_dd {
public:
  cell_incoming_mapping_edge_detection_factor_dd(const INDEX incoming_edge_index, const bool split) : incoming_edge_index_(incoming_edge_index) {}
  cell_incoming_mapping_edge_detection_factor_dd(const INDEX incoming_edge_index) : incoming_edge_index_(incoming_edge_index) {}
  
  template<typename LEFT_FACTOR, typename MSG>
  void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
  {
    assert(msg.size() == l.size());
    for(INDEX i=0; i<msg.size(); ++i) {
      l[i] += msg[i];
    }
  }
  template<typename RIGHT_FACTOR, typename MSG>
  void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
  {
    assert(msg.size() == r.division_distance()-1);
    for(INDEX t=0; t<msg.size(); ++t) {
      r.incoming_transition(incoming_edge_index_, t) += msg[t];
    }
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg -= omega*l;
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    auto diff = r.min_incoming_transition_marginal_diff(incoming_edge_index_);
    msg -= omega*diff;
  } 

  // send messages from detection factor along incoming edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft_deactivated(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    assert(omega <= 1.0 + eps);
    assert(omega > 0.0);
    //std::cout << "send messages to left with omega = " << omega << "\n";

    // check #messages+1 = no incoming edges
    {
      INDEX c=0;
      for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
      assert(c+1 == r.incoming.size());
    }

    const REAL detection_outgoing_cost = r.detection + r.min_outgoing();

    /*
    std::vector<bool> edge_taken(r.incoming.size(), false);
    omega_it = omega_begin;
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = true;
      }
    }
    for(INDEX i=0; i<edge_taken.size()-1; ++i) { edge_taken[i] = true; }


    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    for(INDEX i=0; i<r.incoming.size(); ++i) {
      const REAL val = r.incoming[i];
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    REAL set_to_cost;
    if(smallest_not_taken < smallest_taken) {
      set_to_cost = std::min(detection_outgoing_cost + smallest_not_taken, REAL(0.0));
    } else {
      set_to_cost = std::min({detection_outgoing_cost + second_smallest_taken, detection_outgoing_cost + smallest_not_taken, REAL(0.0)});
    }

    omega_it = omega_begin;
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        const REAL msg = detection_outgoing_cost + r.incoming[ (*msg_it).GetMessageOp().incoming_edge_index_ ] - set_to_cost;
        (*msg_it)[0] -= omega*msg;
      }
    } 
    return;
    */

    const auto smallest_incoming = two_smallest_elements<REAL>(r.incoming.begin(), r.incoming.end()); // do not take into account disappearance cost 

    assert(false); // wrong message computation?
    const REAL set_to_cost = std::min(detection_outgoing_cost + smallest_incoming[1], REAL(0.0));
    //const REAL set_to_cost = std::min(detection_outgoing_cost + smallest_incoming[0], REAL(0.0));
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const REAL msg = detection_outgoing_cost + r.incoming[ (*msg_it).GetMessageOp().incoming_edge_index_ ] - set_to_cost;
      (*msg_it)[1] -= omega*msg;
    } 
  }


  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::vector left_vars,
      RIGHT_FACTOR& r, typename SOLVER::vector detection_vars, 
      typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars,
      typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars ) const
  {
    assert(false);
  }

private:
  const INDEX incoming_edge_index_; 
};

// left is mapping_edge_factor_dd, right is detection_factor_dd
class cell_outgoing_mapping_edge_detection_factor_dd {
public:
  cell_outgoing_mapping_edge_detection_factor_dd(const INDEX outgoing_edge_index, const bool split) : outgoing_edge_index_(outgoing_edge_index) {}
  cell_outgoing_mapping_edge_detection_factor_dd(const INDEX outgoing_edge_index) : outgoing_edge_index_(outgoing_edge_index) {}
  
  template<typename LEFT_FACTOR, typename MSG>
  void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
  {
    assert(msg.size() == l.size());
    for(INDEX i=0; i<msg.size(); ++i) {
      l[i] += msg[i];
    }
  }
  template<typename RIGHT_FACTOR, typename MSG>
  void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
  {
    assert(msg.size() == r.division_distance()-1);
    for(INDEX t=0; t<msg.size(); ++t) {
      r.outgoing_transition(outgoing_edge_index_, t) += msg[t];
    }
    r.outgoing_transition(outgoing_edge_index_, r.outgoing_transition.dim2()-1) += msg[msg.size()-1];
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg -= omega*l;
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    auto diff = r.min_outgoing_transition_marginal_diff(outgoing_edge_index_);
    msg -= omega*diff;
  } 

  // send messages from detection factor along incoming edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft_deactivated(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    assert(omega <= 1.0 + eps);
    assert(omega > 0.0);
    //std::cout << "send messages to left with omega = " << omega << "\n";

    // check #messages+1 = no incoming edges
    {
      INDEX c=0;
      for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
      assert(c+1 == r.outgoing.size());
    }

    const REAL detection_outgoing_cost = r.detection + r.min_outgoing();

    /*
    std::vector<bool> edge_taken(r.incoming.size(), false);
    omega_it = omega_begin;
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = true;
      }
    }
    for(INDEX i=0; i<edge_taken.size()-1; ++i) { edge_taken[i] = true; }


    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    for(INDEX i=0; i<r.incoming.size(); ++i) {
      const REAL val = r.incoming[i];
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    REAL set_to_cost;
    if(smallest_not_taken < smallest_taken) {
      set_to_cost = std::min(detection_outgoing_cost + smallest_not_taken, REAL(0.0));
    } else {
      set_to_cost = std::min({detection_outgoing_cost + second_smallest_taken, detection_outgoing_cost + smallest_not_taken, REAL(0.0)});
    }

    omega_it = omega_begin;
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        const REAL msg = detection_outgoing_cost + r.incoming[ (*msg_it).GetMessageOp().incoming_edge_index_ ] - set_to_cost;
        (*msg_it)[0] -= omega*msg;
      }
    } 
    return;
    */

    const auto smallest_outgoing = two_smallest_elements<REAL>(r.outgoing.begin(), r.outgoing.end()); // do not take into account disappearance cost 

    assert(false); // wrong message computation?
    const REAL set_to_cost = std::min(detection_outgoing_cost + smallest_outgoing[1], REAL(0.0));
    //const REAL set_to_cost = std::min(detection_outgoing_cost + smallest_incoming[0], REAL(0.0));
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const REAL msg = detection_outgoing_cost + r.outgoing[ (*msg_it).GetMessageOp().outgoing_edge_index_ ] - set_to_cost;
      (*msg_it)[1] -= omega*msg;
    } 
  }


  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::vector left_vars,
      RIGHT_FACTOR& r, typename SOLVER::vector detection_vars, 
      typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars,
      typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars ) const
  {
    assert(false);
  } 

private:
  const INDEX outgoing_edge_index_; 
};

// left is division_mapping_edge_factor_dd, right is detection_factor_dd
class cell_incoming_division_edge_detection_factor_dd {
public:
  cell_incoming_division_edge_detection_factor_dd(const INDEX incoming_edge_index, const bool split) : incoming_edge_index_(incoming_edge_index) {}
  cell_incoming_division_edge_detection_factor_dd(const INDEX incoming_edge_index) : incoming_edge_index_(incoming_edge_index) {}
  
  template<typename LEFT_FACTOR>
  void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l.cost += msg;
  }
  template<typename RIGHT_FACTOR>
  void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r.incoming_division[incoming_edge_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg[0] -= omega*l.cost;
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    const REAL diff = r.min_incoming_division_marginal_diff(incoming_edge_index_);
    msg[0] -= omega*diff;
  } 

  // send messages from detection factor along incoming edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft_deactivated(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
  }


  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_var,
      RIGHT_FACTOR& r, typename SOLVER::vector detection_vars, 
      typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars,
      typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars ) const
  {
    assert(false);
  } 
private:
  const INDEX incoming_edge_index_; 
};

// left is division_mapping_edge_factor_dd, right is detection_factor_dd
class cell_outgoing_division_edge_detection_factor_dd {
public:
  cell_outgoing_division_edge_detection_factor_dd(const INDEX outgoing_edge_index, const bool split) : outgoing_edge_index_(outgoing_edge_index) {}
  cell_outgoing_division_edge_detection_factor_dd(const INDEX outgoing_edge_index) : outgoing_edge_index_(outgoing_edge_index) {}
  
  template<typename LEFT_FACTOR>
  void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l.cost += msg;
  }
  template<typename RIGHT_FACTOR>
  void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r.outgoing_division[outgoing_edge_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    msg[0] -= omega*l.cost;
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    const REAL diff = r.min_outgoing_division_marginal_diff(outgoing_edge_index_);
    msg[0] -= omega*diff;
  } 

  // send messages from detection factor along incoming edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft_deactivated(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_var,
      RIGHT_FACTOR& r, typename SOLVER::vector detection_vars, 
      typename SOLVER::vector incoming_division_vars, typename SOLVER::matrix incoming_transition_vars,
      typename SOLVER::matrix outgoing_transition_vars, typename SOLVER::vector outgoing_division_vars ) const
  {
    assert(false);
  } 

private:
  const INDEX outgoing_edge_index_; 
};

class transition_message_dd {
public:
  transition_message_dd(const bool split, const INDEX outgoing_edge_index, const INDEX incoming_edge_index) 
    : 
    outgoing_edge_index_(outgoing_edge_index),
    incoming_edge_index_(incoming_edge_index),
    split_(split)
  {}

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    if(split_) {
      
      msg[0] -= omega*r.min_incoming_division_marginal_diff(incoming_edge_index_);
      return;

    } else {

      auto diff = r.min_incoming_transition_marginal_diff(incoming_edge_index_);
      msg -= omega*diff;
      return;

    }
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  {
    if(split_) {
      
      msg[0] -= omega*l.min_outgoing_division_marginal_diff(outgoing_edge_index_);
      return;

    } else {

      auto diff = l.min_outgoing_transition_marginal_diff(outgoing_edge_index_);
      msg -= omega*diff;
      return;

    }
  }

  // send messages from detection factor along outgoing edges
  template<typename LEFT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToRight(const LEFT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    const auto incoming_transition_min = l.incoming_transition.min2();
    assert(incoming_transition_min.size() == l.division_distance()-1);

    assert(omega > 0.0);
    assert(omega <= 1.0 + eps);

    std::bitset<128> transition_edge_taken(false);
    assert(l.no_outgoing_transition_edges() <= transition_edge_taken.size());
    std::bitset<128> division_edge_taken(false);
    assert(l.no_outgoing_division_edges() <= division_edge_taken.size());

    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      if((*msg_it).GetMessageOp().split_) {
        division_edge_taken[(*msg_it).GetMessageOp().outgoing_edge_index_] = true;
      } else {
        transition_edge_taken[(*msg_it).GetMessageOp().outgoing_edge_index_] = true;
      }
    }

    const INDEX last = l.division_distance()-1;
    REAL smallest_not_taken = 0.0;
    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    const REAL incoming_division_min = l.incoming_division.min();

    auto update_taken = [&smallest_taken, &second_smallest_taken](const REAL x) {
      const REAL min = std::min(smallest_taken, x);
      const REAL max = std::max(smallest_taken, x);
      smallest_taken = min;
      second_smallest_taken = std::min(max, second_smallest_taken); 
    };
    auto update_not_taken = [&smallest_not_taken](const REAL x) { smallest_not_taken = std::min(smallest_not_taken, x); };

    for(INDEX i=0; i<l.no_outgoing_transition_edges(); ++i) {
      if(!transition_edge_taken[i]) {

        update_not_taken( l.detection[0] + incoming_division_min + l.outgoing_transition(i,0) );
        for(INDEX t=1; t<l.division_distance(); ++t) {
          update_not_taken( l.detection[t] + incoming_transition_min[t-1] + l.outgoing_transition(i,t) );
        }

      } else {
        update_taken( l.detection[0] + incoming_division_min + l.outgoing_transition(i,0) );

        for(INDEX t=1; t<l.division_distance()-2; ++t) {
          update_taken( l.detection[t] + incoming_transition_min[t-1] + l.outgoing_transition(i,t) );
        }
        const REAL second_last_val = l.detection[last-1] + incoming_transition_min[last-2] + l.outgoing_transition(i,last-1);
        const REAL last_val = l.detection[last] + incoming_transition_min[last-1] + l.outgoing_transition(i,last);
        update_taken( std::min(second_last_val, last_val) );
      } 
    } 


    for(INDEX i=0; i<l.no_outgoing_division_edges(); ++i) {
      const REAL val = l.detection[last] + incoming_transition_min[last-1] + l.outgoing_division[i];
      if(!division_edge_taken[i]) {
        update_not_taken(val);
      } else {
        update_taken(val);
      }
    }

    const REAL set_to_cost = std::min(smallest_not_taken, second_smallest_taken);
    //const REAL set_to_cost = std::min(smallest_not_taken, smallest_taken);

    assert(std::abs( std::min(smallest_not_taken, smallest_taken) - l.LowerBound()) <= eps);

    vector<REAL> msg(l.division_distance()-1);

    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
        const INDEX edge_index = (*msg_it).GetMessageOp().outgoing_edge_index_;

        // in general, the outgoing division edges need to be read twice, but they may change after first reparametrization. cache them here and revert caching in factors_messages.hxx!
        if((*msg_it).GetMessageOp().split_) {
          (*msg_it)[0] -= 0.5*omega*(l.detection[last] + incoming_transition_min[last-1] + l.outgoing_division[edge_index] - set_to_cost);

        } else { 
          msg[0] = l.detection[0] + incoming_division_min + l.outgoing_transition(edge_index,0);
          for(INDEX t=1; t<l.division_distance()-1; ++t) {
            msg[t] = l.detection[t] + incoming_transition_min[t-1] + l.outgoing_transition(edge_index,t);
          }
          const REAL last_val = l.detection[last] + incoming_transition_min[last-1] + l.outgoing_transition(edge_index,last);
          msg[last-1] = std::min(msg[last-1], last_val);

          for(INDEX i=0; i<msg.size(); ++i) {
            msg[i] -= set_to_cost;
          }

          (*msg_it) -= omega*msg;
        }
      }
  }

  // send messages from detection factor along incoming edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY>
  static void SendMessagesToLeft(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
  {
    const auto outgoing_transition_min = r.outgoing_transition.min2();
    assert(outgoing_transition_min.size() == r.division_distance());

    assert(omega > 0.0);
    assert(omega <= 1.0 + eps);

    std::bitset<128> transition_edge_taken(false);
    assert(r.no_incoming_transition_edges() <= transition_edge_taken.size());
    std::bitset<128> division_edge_taken(false);
    assert(r.no_incoming_division_edges() <= division_edge_taken.size());

    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it) {
      if((*msg_it).GetMessageOp().split_) {
        division_edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = true;
      } else {
        transition_edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = true;
      }
    }

    const INDEX last = r.division_distance()-1; 

    vector<REAL> msg(r.division_distance()-1);
    const REAL outgoing_division_min = r.outgoing_division.min();

    REAL smallest_not_taken = 0.0;
    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();

    auto update_taken = [&smallest_taken, &second_smallest_taken](const REAL x) {
      const REAL min = std::min(smallest_taken, x);
      const REAL max = std::max(smallest_taken, x);
      smallest_taken = min;
      second_smallest_taken = std::min(max, second_smallest_taken); 
    };
    auto update_not_taken = [&smallest_not_taken](const REAL x) { smallest_not_taken = std::min(smallest_not_taken, x); };

    {
      for(INDEX i=0; i<r.no_incoming_division_edges(); ++i) {
        const REAL val = r.detection[0] + outgoing_transition_min[0] + r.incoming_division[i];
        if(division_edge_taken[i]) {
          update_taken(val);
        } else {
          update_not_taken(val);
        }
      } 
    }

    for(INDEX i=0; i<r.no_incoming_transition_edges(); ++i) {
      if(transition_edge_taken[i]) {
        
        for(INDEX t=1; t<r.division_distance(); ++t) {
          update_taken( r.detection[t] + outgoing_transition_min[t] + r.incoming_transition(i,t-1));
        }
        update_taken( r.detection[last] + outgoing_division_min + r.incoming_transition(i,last-1) );

      } else {
        
        for(INDEX t=1; t<r.division_distance(); ++t) {
          update_not_taken( r.detection[t] + outgoing_transition_min[t] + r.incoming_transition(i,t-1) );
        }
        update_not_taken( r.detection[last] + outgoing_division_min + r.incoming_transition(i,last-1));
      }
    }

    assert(std::abs( std::min(smallest_not_taken, smallest_taken) - r.LowerBound()) <= eps);

    const REAL set_to_cost = std::min(smallest_not_taken, second_smallest_taken); 
    //const REAL set_to_cost = std::min(smallest_not_taken, smallest_taken); 

    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
      const INDEX edge_index = (*msg_it).GetMessageOp().incoming_edge_index_;

      if((*msg_it).GetMessageOp().split_) {
        (*msg_it)[0] -= 0.5*omega*(r.detection[0] + r.incoming_division[edge_index] + outgoing_transition_min[0] - set_to_cost);

      } else { 
        for(INDEX t=1; t<r.division_distance()-1; ++t) {
          msg[t-1] = r.detection[t] + outgoing_transition_min[t] + r.incoming_transition(edge_index, t-1);
        }
        msg[last-1] = r.detection[last] + std::min(outgoing_division_min, outgoing_transition_min[last]) + r.incoming_transition(edge_index, last-1);

        for(INDEX i=0; i<msg.size(); ++i) {
          msg[i] -= set_to_cost;
        }

        (*msg_it) -= omega*msg;
      }
    } 
  }

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    assert(split_);
    l.outgoing_division[outgoing_edge_index_] += msg;
    assert(std::isfinite(l.outgoing_division[outgoing_edge_index_])); 
  }

  template<typename G, typename MSG>
  void RepamLeft(G& l, const MSG& msg)
  {
    assert(!split_);
    assert(msg.size() == l.division_distance()-1);
    for(INDEX i=0; i<msg.size(); ++i) {
      l.outgoing_transition(outgoing_edge_index_, i) += msg[i];
      assert(std::isfinite(l.outgoing_transition(outgoing_edge_index_, i)));
    }
    l.outgoing_transition(outgoing_edge_index_, msg.size()) += msg[msg.size()-1];
    assert(std::isfinite(l.outgoing_transition(outgoing_edge_index_, msg.size())));
  }

  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    assert(split_);
    r.incoming_division[incoming_edge_index_] += msg;
    assert(std::isfinite(r.incoming_division[incoming_edge_index_]));
  }

  template<typename G, typename MSG>
  void RepamRight(G& r, const MSG& msg)
  {
    assert(!split_);
    assert(msg.size() == r.division_distance()-1);
    for(INDEX i=0; i<msg.size(); ++i) {
      r.incoming_transition(incoming_edge_index_, i) += msg[i];
      assert(std::isfinite(r.incoming_transition(incoming_edge_index_, i)));
    }
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::vector left_detection_vars, 
      typename SOLVER::vector left_incoming_division_vars, typename SOLVER::matrix left_incoming_transition_vars,
      typename SOLVER::matrix left_outgoing_transition_vars, typename SOLVER::vector left_outgoing_division_vars,
      RIGHT_FACTOR& r, typename SOLVER::vector right_detection_vars, 
      typename SOLVER::vector right_incoming_division_vars, typename SOLVER::matrix right_incoming_transition_vars,
      typename SOLVER::matrix right_outgoing_transition_vars, typename SOLVER::vector right_outgoing_division_vars ) const
  {
    assert(false);
    /*
    if(split_) {
      s.make_equal(l_outgoing_division_literals[outgoing_edge_index_], r_incoming_division_literals[incoming_edge_index_]);

    } else {
      for(INDEX t=0; t<l.division_distance()-2; ++t) {
        s.make_equal( l_outgoing_transition_literals(outgoing_edge_index_, t), r_incoming_transition_literals(incoming_edge_index_, t) );
      }

      const INDEX last = l.division_distance()-1;
      std::array<sat_var,2> last_outgoing_transitions({l_outgoing_transition_literals(outgoing_edge_index_, last-1), l_outgoing_transition_literals(outgoing_edge_index_, last)});
      const auto last_outgoing_literal = s.add_at_most_one_constraint(last_outgoing_transitions.begin(), last_outgoing_transitions.end());
      s.make_equal(last_outgoing_literal, r_incoming_transition_literals(incoming_edge_index_, last-1));
    }
    */
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    constexpr INDEX no_edge_taken = detection_factor_dd::no_edge_taken;
    const auto& lp = l.primal();
    const auto& rp = r.primal();
    return true;

    if(split_) {
      const bool left_active = lp.outgoing_division && (lp.outgoing_edge == outgoing_edge_index_);
      const bool right_active = rp.division == 0 && (rp.incoming_edge == incoming_edge_index_);
      if(left_active != right_active) {
        std::cout << "division edge inconsistency.\n";
      }
      return left_active == right_active;
    } else {
      const bool left_active = (lp.outgoing_division == false) && (lp.outgoing_edge == outgoing_edge_index_);
      const bool right_active = (rp.division > 0) && (rp.incoming_edge == incoming_edge_index_); 
      if(left_active && right_active) {
        if( std::min(lp.division+1, l.division_distance()-1) != rp.division) {
          std::cout << "transition edge inconsistency: both edges on, but division distance not respected.\n";
        }
        return std::min(lp.division+1, l.division_distance()-1) == rp.division;
      } else {
        if(left_active != right_active) {
          std::cout << "transition edge inconsistency in non-last edge.\n";
        }
        return left_active == right_active;
      }
    }
  }

private:
  const INDEX outgoing_edge_index_;
  const INDEX incoming_edge_index_;
  const bool split_;
};

} // end namespace LP_MP

#endif // LP_MP_DETECTION_FACTOR_HXX
