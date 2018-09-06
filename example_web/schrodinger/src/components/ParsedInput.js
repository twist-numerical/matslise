import React, { Component } from 'react';
import DF from '../lib/Differentiable';

class ParsedInput extends Component {
  static defaultProps = {
    onParsed: () => {},
    parseOnMount: true,
  }
  state = {
    value: "",
    parsed: null,
  }
  input = null;

  constructor(props) {
    super(props);

    this.state.value = props.value;
  }

  componentDidMount() {
    if(this.props.parseOnMount)
      this.props.onParsed(this.onInput(this.state.value));
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if(prevProps.value !== this.props.value) {
      this.onInput(this.props.value); 
    }
  }

  render() {
    return (
      <div className="input-group mb-3">
      <div className="form-control" contentEditable={true}
      style={{'textAlign': 'left'}}
      onInput={(e) => this.onInput()}
      ref={div => this.input = div} />
      <div className="input-group-append">
      <button
      onClick={() => this.props.onParsed(this.state.parsed)}
      className="btn btn-outline-secondary">View</button>
      </div>
      </div>
      );
  }

  getCursorOffset(node) {
    if(!node)
      return 0;

    if (window.getSelection !== undefined) {
      const selection = window.getSelection();
      if(selection.type !== "None") {
        const range = selection.getRangeAt(0);
        const preCaretRange = range.cloneRange();
        preCaretRange.selectNodeContents(node);
        preCaretRange.setEnd(range.endContainer, range.endOffset);
        return preCaretRange.toString().length;
      }
    } else if (document.selection !== undefined
      && document.selection.type !== "Control"
      && document.selection.type !== "None") {
      const textRange = document.selection.createRange();
      const preCaretTextRange = document.body.createTextRange();
      preCaretTextRange.moveToElementText(node);
      preCaretTextRange.setEndPoint("EndToEnd", textRange);
      return preCaretTextRange.text.length;
    }

    return 0;
  }

  setCursor(node, position){
    if(!node)
      return;

    if(document.createRange) {
      const range = document.createRange();
      range.selectNodeContents(node);
      range.setStart(node, position);
      range.setEnd(node, position);
      const selection = window.getSelection();
      selection.removeAllRanges();
      selection.addRange(range);
    } else if(node.createTextRange) {
      const textRange = node.createTextRange();
      textRange.collapse(true);
      textRange.moveEnd(position);
      textRange.moveStart(position);
      textRange.select();
    } else if(node.setSelectionRange) {
      node.setSelectionRange(position,position);
    }
  }

  onInput(value) {
    if(value === undefined)
      value = this.input.textContent;
    const {parsed, todo, done} = this.parse(value);
    const cursor = this.getCursorOffset(this.input);
    this.setState({parsed, value});
    
    const doneSpan = document.createElement('span');
    doneSpan.textContent = done;
    
    const todoSpan = document.createElement('span');
    todoSpan.style.color = 'red';
    todoSpan.textContent = todo;
    
    while (this.input.firstChild)
      this.input.removeChild(this.input.firstChild);
    this.input.appendChild(doneSpan);
    this.input.appendChild(todoSpan);

    if(cursor <= done.length)
      this.setCursor(doneSpan.firstChild, cursor);
    else
      this.setCursor(todoSpan.firstChild, cursor-done.length);

    return parsed;
  }

  parse(value) {
    let parsed = null;
    let done = value;
    let todo = "";
    try {
      parsed = DF.parse(value);
    } catch (e) {
      if(e.done === undefined || e.todo === undefined)
        throw e;
      done = e.done;
      todo = e.todo;
    }
    return {
      parsed: parsed,
      done: done,
      todo: todo,
    };
  }
}

export default ParsedInput;
