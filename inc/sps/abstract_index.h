#pragma once


namespace sps
{

/**
 * @brief Abstract Superclass for all Indices.
 *
 * Not usefull on it's own.
 */
class AbstractIndex
{
  public:
    /**
     * @brief Does nothing.
     * required to make AbstractIndex downcastable to Index<type_defs> in python by making it abstract
     */
    virtual void dummy( )
    {}
};


} // namespace sps
