---
name: deep-planning
description: Deep, context-first interrogation of a plan or design before any implementation - you carry the context load, the user answers in one line, and the output is a spec including its decision record. Invoke ONLY when the user names it - "deep planning", "use deep planning", "deep plan this", "用 deep planning", "深入规划", "grill me on this", or /deep-planning. Do NOT invoke it for ordinary planning, design, or "let's build X" requests where the user has not named it; superpowers:brainstorming handles those.
---

# Deep Planning

## Why this exists

In a design conversation the expensive failure is not the user typing a long answer. It is you
locking in a wrong assumption because the question was never properly framed.

Ordinary brainstorming skills shorten the question to save the user effort — "one question at a
time", "multiple choice preferred". That inverts the wrong thing. A one-line question forces the
user to reconstruct, unaided, what you already know, what you have assumed, and what the options
even are. The effort didn't disappear; it moved onto them, in the form they are worst placed to
supply.

This skill puts the context load back on you.

---

## First principle: asymmetric effort

**The context is your job. The decision is theirs.**

Every question you ask MUST carry, in the question itself:

- the current state, and the facts you have already looked up
- the options, 2–4 of them, each with its trade-off
- your recommended answer, with the reasoning
- where your recommendation could be wrong

Every answer they give MAY be one line. **"B" is a complete and valid answer.** Never ask the user
to justify their choice. Never respond to a short answer by re-asking the same question in other
words. If a short answer leaves something genuinely unresolved, name the specific consequence you
are unsure about — do not hand the question back.

---

## Standing rules

These apply for the entire session. This file is read once and not re-read on later turns, so
treat what follows as continuously-in-force behaviour, not a checklist you complete and discard.

1. **No question before the context brief.** (Phase 0.)
2. **Facts you can look up, you look up.** Only *decisions* go to the user.
3. **Every check-list item gets one of three fates:** asked, answered by you from the repo, or
   declared not-applicable with a reason. Never silently skipped.
4. **Never close a question silently.** Say the call and its basis out loud, every time.
5. **Maintain the running log, and re-read it before composing each question.** Your grip on
   these rules decays over a long session; the log is the anchor.

---

## Phase 0 — Context sweep. Ask nothing yet.

Before the first question, go and find out. Read the repo structure, recent specs and plans, the
labbook, recent commits, the files the request touches. Work through the check list below and
resolve every item you can resolve without the user.

Then produce a **context brief** and nothing else:

> **What I found** — the facts that bear on this, with paths.
> **What I'm assuming** — assumptions you are proceeding on unless corrected.
> **What's actually undecided** — the short list that genuinely needs them.

**Hard gate: no questions until the context brief has been delivered.** A brief that says "I
couldn't find much" is fine; skipping the brief is not.

### Scale gate

Judge the size of the work first. If it is small or mechanical, say so plainly — *"this doesn't
need the full treatment; here's what I'd do, say go"* — and let them wave you through. Grinding a
two-file change through a full interrogation is its own failure.

---

## Phase 1 — Lay out the decision tree

List every open decision, ordered by dependency: the ones that constrain other decisions come
first. Show the user **the whole list at once**, before asking question one.

They should know how many questions are coming and what they are. Being fed one question at a
time with no idea how far the road runs is exactly the experience this skill exists to replace.

If answering an early question collapses or adds later ones, say so when it happens and show the
revised list.

---

## Phase 2 — The question loop

One question per message. Each one in this shape:

```
【Context】   why this question arises now; the facts you found that bear on it
【Question】  one, specific
【Options】   2–4, each with its trade-off
【My recommendation】  + the reasoning
【Where I could be wrong】  the condition under which your recommendation is the wrong call
      ↓  they answer (one line is fine)
【Drill-down】  only if the decision changes the implementation
      ↓
【Closing】   "Closing this one — <basis>. Say so if you disagree."
```

### When to drill down

The test: **would answering A rather than B produce different code, or a different experiment?**

- **No** — close it in one round. It is a preference, not a fork.
- **Yes** — drill into the specific consequence: what it forces downstream, what it rules out,
  whether it contradicts a decision already made.

This judgement is contextual and you will sometimes get it wrong. That is acceptable. What is not
acceptable is making it silently.

### Closing out loud

Before moving on, state the call and its basis in one sentence:

> "Closing this — A and B both end up running the same experiment, so it doesn't reach
> implementation. Say so if you disagree."

The user then overrides with one line, or says nothing and you proceed. Soft judgement, visible
reasoning, cheap override. **A silent "that's enough, next question" is the failure mode this
whole skill is built to prevent** — it is how thorough questioning decays into thin questioning
without anyone noticing.

### Decisions you make on their behalf

For low-stakes points, decide yourself rather than spending a question. Say which way you went
and why, in one line, and carry on: *"Called this myself: one check list, not two — say if you'd
rather split it."*

---

## Phase 3 — Write it up

One document. Default location `docs/superpowers/specs/YYYY-MM-DD-<topic>-design.md`, unless the
project's own conventions say otherwise.

Normal spec content — motivation, scope, architecture, testing — **plus** a decision record:

| Decision | What was decided | Why | Rejected, and why | What would overturn it |
| --- | --- | --- | --- | --- |

The last column is the one that earns its keep. It marks a decision's shelf life, so the next
session can tell a live constraint from a stale one instead of re-litigating it or silently
violating it.

Then hand off: **`superpowers:writing-plans`**.

---

## The check list

One list. Items apply as relevant — do not sort the work into "code" or "experiment" first and
then ask only half. Each item is asked, answered by you, or declared not-applicable with a reason.

**Intent**
- What is the hypothesis, and what result would falsify it?
- What is the baseline or control? (An experiment without one does not start.)
- What is the success criterion — fixed **before** the run, not chosen after seeing the numbers?
- What happens on a positive, a negative, and an ambiguous result? Does a negative result still
  get written up?

**Validity**
- What are the known confounds? Which are controlled, which are knowingly accepted?
- What existing code, splits, or checkpoints get reused, and does reusing them leak anything?

**Cost**
- What is the compute budget, and what gets cut first if it overruns?

**Build**
- Which files are created or modified, and where do the boundaries fall?
- Does this break an existing interface?
- How is it tested — what specifically is asserted?
- What counts as done?

---

## Red flags

Each of these thoughts means you are about to under-question. Stop.

| Thought | Reality |
| --- | --- |
| "This is obvious, I'll just assume it" | Assumptions are the failure this skill exists to catch. Ask, or state the assumption in the brief. |
| "I'll ask, it's faster than looking it up" | Facts are your job. Look it up. Only decisions go to them. |
| "They gave a short answer, they must want to move on" | Short answers are *permitted by design*. Length says nothing about whether the question is closed. |
| "I'll skip the brief, I already know this repo" | The brief is how they check whether you know it. Write it. |
| "That check-list item doesn't apply here" | Fine — say so, with the reason. Silence is not a fate. |
| "Enough on this one, next question" | Not unless you said it out loud and they had the chance to object. |
| "I'll ask the remaining three together to save time" | One per message. Bundled questions get one answer and two silent assumptions. |
