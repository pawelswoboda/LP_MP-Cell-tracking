#ifndef LP_MP_DETECTION_FACTOR_FINE_FACTORS_HXX
#define LP_MP_DETECTION_FACTOR_FINE_FACTORS_HXX

namespace LP_MP {

  // used for incoming and outgoing edges
  class flow_conservation_factor {
    public:
      flow_conservation_factor(const INDEX size) 
        : edges(size,0.0),
        detection(0.0)
    {}

      REAL min_detection_cost() const { return edges.min() + detection; }
      REAL LowerBound() const { return std::min(0.0, min_detection_cost()); }
      REAL EvaluatePrimal() const 
      {
        if(primal != std::numeric_limits<INDEX>::max()) {
          assert(primal < edges.size());
          return detection + edges[primal];
        } else {
          return 0.0;
        }
      }

      REAL min_marginal(const INDEX i) const
      {
        assert(i<edges.size());
        const REAL cost_active = detection + edges[i];
        const REAL orig_edge_val = edges[i];
        edges[i] = std::numeric_limits<REAL>::infinity();
        const REAL cost_inactive = std::min(0.0, detection + edges.min());
        edges[i] = orig_edge_val;
        return cost_active - cost_inactive;
      }

  void set_incoming_transition_cost(const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(edges[edge_index] == 0.0); 
    edges[edge_index] = cost; 
  } 
  void set_outgoing_transition_cost(const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(edges[edge_index] == 0.0);
    edges[edge_index] = cost;
  }
  void set_incoming_division_cost(const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(edges[no_incoming_transition_edges + edge_index] == 0.0);
    edges[no_incoming_transition_edges + edge_index] = cost; 
  } 
  void set_outgoing_division_cost(const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, const INDEX edge_index, const REAL cost) 
  { 
    assert(edges[no_outgoing_transition_edges + edge_index] == 0.0);
    edges[no_outgoing_transition_edges + edge_index] = cost; 
  }


  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( detection, edges ); }
  auto export_variables() { return std::tie( detection, edges ); }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::variable detection_var, typename SOLVER::vector edges_var) const
  {
    assert(false);
  }
  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::variable detection_var, typename SOLVER::vector edges_var) const
  {
    assert(false);
  }
  
  void init_primal() { primal = std::numeric_limits<INDEX>::max(); }
  bool detection_active() const { return detection < edges.size(); }
  


    mutable vector<REAL> edges;
    REAL detection;
    INDEX primal; // if std::numeric_limits<INDEX>::max, then not active
  };

  // message between two flow_conservation factor belong to subsequent time steps
  class flow_conservation_message {
    public:
    flow_conservation_message(const INDEX outgoing_edge_index, const INDEX incoming_edge_index, const bool split) 
      : 
        outgoing_edge_index_(outgoing_edge_index),
        incoming_edge_index_(incoming_edge_index),
        split_(split)
    {
      assert(outgoing_edge_index <= std::numeric_limits<SHORT_INDEX>::max());
      assert(incoming_edge_index <= std::numeric_limits<SHORT_INDEX>::max());
    }

    template<typename RIGHT_FACTOR, typename G2>
    void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
    {
      msg[0] -= omega*r.min_marginal(incoming_edge_index_);
    }

    template<typename LEFT_FACTOR, typename G2>
    void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
    {
      msg[0] -= omega*l.min_marginal(outgoing_edge_index_);
    }

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l.edges[outgoing_edge_index_] += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r.edges[incoming_edge_index_] += msg;
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_var, typename SOLVER::vector left_edges,
      RIGHT_FACTOR& r,  typename SOLVER::variable right_var, typename SOLVER::vector right_edges) const
  {
    assert(false);
  }

    private:
      const SHORT_INDEX outgoing_edge_index_;
      const SHORT_INDEX incoming_edge_index_;
      bool split_;
  };

  // message between detection variable
  class detection_message {
    public:

    template<typename RIGHT_FACTOR, typename G2>
    void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
    {
      msg[0] -= omega*r.min_detection_cost();
    }

    template<typename LEFT_FACTOR, typename G2>
    void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
    {
      msg[0] -= omega*l.min_detection_cost();
    }

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
    r.detection += msg;
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, 
      LEFT_FACTOR& l, typename SOLVER::variable left_var, typename SOLVER::vector left_edges,
      RIGHT_FACTOR& r,  typename SOLVER::variable right_var, typename SOLVER::vector right_edges) const
  {
    assert(false);
  }
  };
} // namespace LP_MP

#endif // LP_MP_DETECTION_FACTOR_FINE_FACTORS_HXX
