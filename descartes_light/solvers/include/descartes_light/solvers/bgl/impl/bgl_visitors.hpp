#ifndef DESCARTES_LIGHT_BOOST_VISITORS
#define DESCARTES_LIGHT_BOOST_VISITORS

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <descartes_light/solvers/bgl/boost_ladder_types.h>
#include <descartes_light/core/solver.h>
#include <descartes_light/core/edge_evaluator.h>
#include <descartes_light/solvers/bgl/boost_ladder_types.h>

//heres a terrible idea: pass the graph by reference into the vsitor so that it may edit as it traverses
namespace descartes_light
{
template <typename FloatType>
struct AddAllVisitor : boost::default_dijkstra_visitor
{
  AddAllVisitor(std::vector<typename EdgeEvaluator<FloatType>::ConstPtr>& edge_eval,
          std::map<VertexDesc<FloatType>, VertexDesc<FloatType>>& predecessor_map,
          std::vector<std::vector<VertexDesc<FloatType>>>& ladder_rungs,
          bglgraph<FloatType>& eg)
          : eval(edge_eval),
            predecessors(predecessor_map),
            rungs(ladder_rungs),
            mutable_graph(eg)
       {
       }

  // must have copy constructor
  //AddAllVisitor() : AddAllVisitor

  void examine_vertex(VertexDesc<FloatType> u, bglgraph<FloatType> g)
  {
    int out_deg = boost::out_degree(u, g);
    // return if the vertex has any out edges
    if (out_deg == 0)
    {
      // check rung
      VertexDesc<FloatType> curr = u;
      // todo: find a better way to select rungs
      std::size_t next_rung = 0; // assuming that a dummy start node is still used -> This does not work bc predecessors preallocates based on vertices
      for (VertexDesc<FloatType> prev = predecessors.at(curr); prev != curr; curr = prev, prev = predecessors.at(curr))
      {
        ++next_rung;
      }
      if (next_rung < rungs.size()) //compiler will optimize more if this branching could be avoided
      {
        FloatType cost;
        for (std::size_t s = 0; s < rungs[next_rung].size(); ++s)
        {
          std::pair<bool, FloatType> results =
              eval[static_cast<size_t>(next_rung-1)]->evaluate(*g[u].state, *g[rungs[next_rung][s]].state);
          if (results.first)
          {
            cost = results.second + g[rungs[next_rung][s]].cost;
            if (next_rung == 1) //this if can probably be captured in the dummy node edges
              cost += g[u].cost;
            VertexDesc<FloatType> sap = rungs[next_rung][s];
            boost::add_edge(u, sap, cost, mutable_graph);
          }
        }
      }
    }
    return;
  }

private:
  std::vector<typename EdgeEvaluator<FloatType>::ConstPtr> eval{ nullptr };
  std::map<VertexDesc<FloatType>, VertexDesc<FloatType>>& predecessors;
  std::vector<std::vector<VertexDesc<FloatType>>>& rungs;
  bglgraph<FloatType>& mutable_graph{ nullptr };
};
} //descartes_light
#endif
