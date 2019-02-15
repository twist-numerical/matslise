import React, { Component } from 'react';
import DF from '../lib/Differentiable';

const htmlEscape = (s="") => {
  return s.replace(/&/g, '&amp;')
  .replace(/>/g, '&gt;')
  .replace(/</g, '&lt;')
  .replace(/"/g, '&quot;')
  .replace(/'/g, '&apos;');
}

class ParsedInput extends Component {

  static defaultProps = {
    onParsed: () => {},
    onEnter: (pi) => {
      let e = pi.input;
      while(e && e.tagName !== "FORM")
        e = e.parentElement;
      if(e) {
        let submit = document.createElement("input");
        submit.type = "submit";
        submit.style.display = "none";
        e.appendChild(submit);
        submit.click();
        e.removeChild(submit);
      }
    },
    parsed: null,
    variables: [],
  }
  state = {
    parsed: this.props.parsed,
    done: this.props.parsed ? this.props.parsed.value : "", 
    todo: "",
  }
  input = null;

  isValid() {
    return this.state.parsed !== null;
  }

  render() {
    const {parsed, onParsed, variables, onEnter, ...restProps} = this.props;

    return (
      <div {...restProps} contentEditable={true}
      onInput={e => this.onInput(e.target.textContent.replace(/[\n\r]/g, " "))}
      ref={div => this.input = div}
      dangerouslySetInnerHTML={{__html:
        '<span>'+htmlEscape(this.state.done)+'</span><span style="color: red;">'+htmlEscape(this.state.todo)+'</span>'}}
      onKeyDown={(e) => {
        if (e.keyCode === 13) {
          this.props.onEnter(this);
          e.preventDefault();
        }
      }}
      />
      );
  }

  componentDidUpdate() {
    if(!this.input.firstChild)
      return;
    if(this.cursorPosition >= 0) {
      if (this.cursorPosition <= this.state.done.length)
        this.setCursor(this.input.children[0].firstChild, this.cursorPosition);
      else
        this.setCursor(this.input.children[1].firstChild, this.cursorPosition-this.state.done.length);
    }
    this.cursorPosition = -1;
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
    this.cursorPosition = this.getCursorOffset(this.input);

    const {parsed, todo, done} = this.parse(value);

    if (parsed !== null)
      this.props.onParsed(parsed);

    this.setState({parsed, todo, done});
  }

  parse(value) {
    let parsed = null;
    let done = value;
    let todo = "";
    try {
      parsed = DF.parse(value, this.props.variables);
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
