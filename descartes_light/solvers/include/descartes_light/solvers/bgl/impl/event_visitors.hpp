#ifndef DESCARTES_LIGHT_SOLVERS_BGL_IMPL_EVENT_VISITORS_HPP
#define DESCARTES_LIGHT_SOLVERS_BGL_IMPL_EVENT_VISITORS_HPP

#include <descartes_light/solvers/bgl/boost_graph_types.h>

#include <descartes_light/descartes_macros.h>
DESCARTES_IGNORE_WARNINGS_PUSH
#include <boost/graph/visitors.hpp>
DESCARTES_IGNORE_WARNINGS_POP

namespace descartes_light
{
/**
 * @brief Event visitor that terminates the search when a vertex in the last rung of the graph is encountered
 * @details Throws the vertex descriptor that is the termination of the path once a vertex in the last rung of
 * the graph is encountered
 */
template <typename FloatType>
struct early_terminator : public boost::base_visitor<early_terminator<FloatType>>
{
  /** @brief Event filter typedef defining the events for which this visitor can be used */
  typedef boost::on_examine_vertex event_filter;

  early_terminator(long last_rung_idx) : last_rung_idx_(last_rung_idx) {}

  void operator()(VertexDesc<FloatType> u, const BGLGraph<FloatType>& g)
  {
    if (g[u].rung_idx == last_rung_idx_)
      throw u;
  }

  const long last_rung_idx_;
};


/**
 * @brief Event visitor that adds all edges to each node discovered on a depth first search
 * @details Throws the vertex descriptor that is the termination of the path once a vertex in the last rung of
 * the graph is encountered
 */
template <typename FloatType>
struct add_all_dfs : public boost::base_visitor<add_all_dfs<FloatType>>
{
  /** @brief Event filter typedef defining the events for which this visitor can be used */
  typedef boost::on_discover_vertex event_filter;

  add_all_dfs(std::vector<typename EdgeEvaluator<FloatType>::ConstPtr>& edge_eval,
               std::vector<std::vector<VertexDesc<FloatType>>>& ladder_rungs)
               : eval_(edge_eval),
                 ladder_rungs_(ladder_rungs)
  {}

  void operator()(VertexDesc<FloatType> u, const BGLGraph<FloatType>& g)
    {
      auto out_deg = boost::out_degree(u, g);
      // return if the vertex has any out edges
      if (out_deg == 0)
      {
        std::size_t current_rung = static_cast<std::size_t>(g[u].rung_idx);
        std::size_t next_rung = static_cast<std::size_t>(current_rung + 1);
        if (next_rung < ladder_rungs_.size())
        {
          FloatType cost;
          for (std::size_t s = 0; s < ladder_rungs_[next_rung].size(); ++s)
          {
            std::pair<bool, FloatType> results =
                eval_[static_cast<size_t>(current_rung)]->evaluate(*g[u].sample.state, *g[ladder_rungs_[next_rung][s]].sample.state);
            if (results.first)
            {
              cost = results.second + g[ladder_rungs_[next_rung][s]].sample.cost;
              if (current_rung == 0)
                cost += g[u].sample.cost;
              VertexDesc<FloatType> sap = ladder_rungs_[next_rung][s];
              BGLGraph<FloatType>* mutable_graph_ = const_cast<BGLGraph<FloatType>*>(&g);
              boost::add_edge(u, sap, cost, *mutable_graph_);
            }
          }
        }
        // Depth First Component
        else
        {
          throw u;
        }
      }
      return;
    }

  std::vector<typename EdgeEvaluator<FloatType>::ConstPtr> eval_{ nullptr };
  std::vector<std::vector<VertexDesc<FloatType>>>& ladder_rungs_;
};

/**
 * @brief Event visitor for updating vertex cost
 */
template <typename FloatType>
struct cost_recorder : public boost::base_visitor<cost_recorder<FloatType>>
{
  /** @brief Event filter typedef defining the events for which this visitor can be used */
  typedef boost::on_tree_edge event_filter;

  cost_recorder()
  {}

  void operator()(EdgeDesc<FloatType> e, const BGLGraph<FloatType>& g)
  {
    VertexDesc<FloatType> target = boost::target(e, g);
    VertexDesc<FloatType> source = boost::source(e, g);
    BGLGraph<FloatType>* mutable_graph_ = const_cast<BGLGraph<FloatType>*>(&g);

    mutable_graph_->operator[](target).distance = g[source].distance + g[target].sample.cost;
  }
};
}  // namespace descartes_light

#endif  // DESCARTES_LIGHT_SOLVERS_BGL_IMPL_EVENT_VISITORS_HPP